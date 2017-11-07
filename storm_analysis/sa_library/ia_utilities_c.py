#!/usr/bin/env python
"""
Some image analysis utility functions, such as finding local
maxima and identifying peaks that are near each other.

Some of the lifting is done by ia_utilities.c.

Hazen 10/17
"""
import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import scipy
import scipy.spatial

import storm_analysis.sa_library.loadclib as loadclib

util = loadclib.loadCLibrary("storm_analysis.sa_library", "ia_utilities")


# Module constants, these must match the corresponding constants in multifit.h.
RUNNING = 0
CONVERGED = 1
ERROR = 2

# C flmData structure.
class flmData(ctypes.Structure):
    _fields_ = [('margin', ctypes.c_int),
                ('npeaks', ctypes.c_int),
                ('radius', ctypes.c_int),
                ('zrange', ctypes.c_int),

                ('xsize', ctypes.c_int),
                ('ysize', ctypes.c_int),
                ('zsize', ctypes.c_int),

                ('threshold', ctypes.c_double),
                
                ('zvalues', ctypes.c_void_p),

                ('taken', ctypes.c_void_p),
                ('images', ctypes.c_void_p)]

# C interface definition.
util.calcMaxPeaks.argtypes = [ctypes.POINTER(flmData)]
util.calcMaxPeaks.restype = ctypes.c_int

util.findLocalMaxima.argtypes = [ctypes.POINTER(flmData),
                                 ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.float64)]


class MaximaFinder(object):
    """
    For finding local maxima in an image or an image stack.
    """
    def __init__(self, margin = None, n_duplicates = 1, radius = None, threshold = None, z_range = None, z_values = None, **kwds):
        """
        margin - Margin around the edge to avoid in the search.
        n_duplicates - How many times we can return the same location (in multiple calls) before
                       giving up.
        radius - Radius in pixels over which the current pixel needs to be a maxima.
        threshold - Minimum value to be considered as a maxima.
        z_range - Number of z planes over which the current pixel needs to be a maxima.
        z_values - Z values to fill in the z array, there should be one of these for each plane.
        """
        super(MaximaFinder, self).__init__(**kwds)

        self.n_duplicates = n_duplicates
        self.n_planes = len(z_values)
        self.taken = None

        self.p_data = ctypes.c_void_p * self.n_planes
        self.c_taken = self.p_data()

        # Use full range if z_range is not specified.
        if z_range is None:
            z_range = self.n_planes
        
        # Create flmData structure that we'll pass to the C functions
        # that will do all the heavy lifting.
        #
        self.c_zvalues = numpy.ascontiguousarray(numpy.array(z_values), dtype = numpy.float64)
        self.flm_data = flmData(margin = int(margin),
                                radius = int(radius),
                                zrange = int(z_range),
                                threshold = float(threshold),
                                zsize = self.n_planes,
                                zvalues = self.c_zvalues.ctypes.data)

    def findMaxima(self, images):
        """
        Find the (local) maxima in a list of one or more images, assumed to be in the
        same order as z_values.

        Note: This will destructively modify the images, make a copy if you don't want
              them changed.
        """
        assert (len(images) == self.n_planes), "Number of planes does not match number of Z planes."

        # Create 'taken' arrays if they don't already exist. We'll use these to keep
        # track of whether or not we have previously returned a particular pixel.
        #
        if self.taken is None:
            self.taken = []
            for i in range(self.n_planes):
                taken = numpy.ones(images[0].shape, dtype = numpy.int32) - self.n_duplicates
                taken = numpy.ascontiguousarray(taken)
                self.c_taken[i] = taken.ctypes.data
                self.taken.append(taken)
            self.flm_data.taken = ctypes.c_void_p(ctypes.addressof(self.c_taken))

            # Also set the image size.
            self.flm_data.xsize = int(images[0].shape[1])
            self.flm_data.ysize = int(images[0].shape[0])

        # Verify images are the correct size and C contiguous.
        for image in images:
            assert (image.shape[0] == self.taken[0].shape[0]), "Unexpected image x size!"
            assert (image.shape[1] == self.taken[0].shape[1]), "Unexpected image y size!"
            assert (image.flags['C_CONTIGUOUS']), "Image is not C contiguous!"

        # Create pointer array to the images.
        c_images = self.p_data()
        for i in range(self.n_planes):
            c_images[i] = images[i].ctypes.data

        # Set flm_data images pointer.
        self.flm_data.images = ctypes.c_void_p(ctypes.addressof(c_images))

        # Figure out the maximum possible number of peaks.
        max_npeaks = util.calcMaxPeaks(ctypes.byref(self.flm_data))

        # Allocate storage for x,y,z locations
        c_x = numpy.ascontiguousarray(numpy.zeros(max_npeaks, dtype = numpy.float64))
        c_y = numpy.ascontiguousarray(numpy.zeros(max_npeaks, dtype = numpy.float64))
        c_z = numpy.ascontiguousarray(numpy.zeros(max_npeaks, dtype = numpy.float64))

        self.flm_data.npeaks = max_npeaks

        # Get peak locations.
        util.findLocalMaxima(ctypes.byref(self.flm_data), c_z, c_y, c_x)

        np = self.flm_data.npeaks
        return [c_x[:np], c_y[:np], c_z[:np]]

    def resetTaken(self):
        """
        Restore the taken arrays to their original values.
        """
#        # This 
#        if self.taken is None:
#            return
        
        for i, taken in enumerate(self.taken):
            taken[:,:] = 1 - self.n_duplicates

            # Check that the above did not move the array.
            assert (taken.ctypes.data == self.c_taken[i])

def runningIfHasNeighbors(status, c_x, c_y, n_x, n_y, radius):
    """
    Update status based on proximity of new peaks (n_x, n_y) to current peaks (c_x, c_y).

    This works the simplest way by making a KD tree from the new peaks then comparing
    the old peaks against this tree. However this might not be the fastest way given
    that there will likely be a lot more current peaks then new peaks.
    """
    # Make kdtree from new peaks.
    kd = scipy.spatial.KDTree(numpy.stack((n_x, n_y), axis = -1))

    # Test (converged) old peaks for proximity to new peaks.
    for i in range(status.size):
        if (status[i] == CONVERGED):
            [dist, index] = kd.query(numpy.array([[c_x[i], c_y[i]]]), distance_upper_bound = radius)
            if (dist[0] != numpy.inf):
                status[i] = RUNNING

    return status
        

#
# The MIT License
#
# Copyright (c) 2017 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
