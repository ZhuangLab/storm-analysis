#!/usr/bin/env python
"""
Utility classes functions that are used for drift correction.

Hazen 02/17
"""

import numpy
import scipy

import storm_analysis.sa_library.grid_c as gridC
import storm_analysis.sa_library.sa_h5py as saH5Py


class SAH5DriftCorrection(saH5Py.SAH5Py):
    """
    A sub-class of SAH5Py designed for use in drift correction.
    """
    def __init__(self, scale = None, z_bins = 1, **kwds):
        super(SAH5DriftCorrection, self).__init__(**kwds)
        
        self.dx = 0.0
        self.dy = 0.0
        self.dz = 0.0
        self.fmin = None
        self.fmax = None
        self.im_shape_2D = (self.hdf5.attrs['movie_x']*scale,
                            self.hdf5.attrs['movie_y']*scale)
        self.im_shape_3D = (self.hdf5.attrs['movie_x']*scale,
                            self.hdf5.attrs['movie_y']*scale,
                            z_bins)
        self.scale = scale
        self.z_bins = z_bins
        
    def grid2D(self, drift_corrected = False):
        image = numpy.zeros(self.im_shape_2D, dtype = numpy.int32)
        for locs in self.locsInFrameRangeIterator(self.fmin, self.fmax, ["x", "y"]):
            if drift_corrected:
                locs["x"] += self.dx
                locs["y"] += self.dy
            i_x = numpy.floor(locs["x"]*self.scale).astype(numpy.int32)
            i_y = numpy.floor(locs["y"]*self.scale).astype(numpy.int32)
            gridC.grid2D(i_x, i_y, image)
        return image

    def grid3D(self, z_min, z_max, drift_corrected = False):
        z_scale = float(self.z_bins)/(z_max - z_min)
        image = numpy.zeros(self.im_shape_3D, dtype = numpy.int32)
        for locs in self.locsInFrameRangeIterator(self.fmin, self.fmax, ["x", "y", "z"]):

            # Create z value filter.
            #
            # We filter here rather than just relying on gridC.grid3D as application
            # of z drift correction could move out of z range peaks into the acceptable
            # range.
            #
            mask = (locs["z"] > z_min) & (locs["z"] < z_max)
            if (numpy.count_nonzero(mask) == 0):
                continue

            # Remove localizations that are out of range.
            locs["x"] = locs["x"][mask]
            locs["y"] = locs["y"][mask]
            locs["z"] = locs["z"][mask]
            
            # Apply drift correction if requested.
            if drift_corrected:
                locs["x"] += self.dx
                locs["y"] += self.dy
                locs["z"] += self.dz
                        
            # Add to image.
            i_x = numpy.floor(locs["x"]*self.scale).astype(numpy.int32)
            i_y = numpy.floor(locs["y"]*self.scale).astype(numpy.int32)
            i_z = numpy.floor((locs["z"] - z_min)*z_scale).astype(numpy.int32)
            gridC.grid3D(i_x, i_y, i_z, image)
        return image

    def locsInFrameRangeIterator(self, start, stop, fields):
        for i in range(start, stop):
            locs = self.getLocalizationsInFrame(i,
                                                drift_corrected = False,
                                                fields = fields)

            # Skip empty frames.
            if not bool(locs):
                continue
            
            yield locs

    def saveDriftData(self, all_dx, all_dy, all_dz):
        """
        Store drift correction data in the HDF5 file. The all_** arrays
        contain the drift corrections for every frame in the movie in
        units of pixels (X,Y) or microns (Z).
        """
        assert(len(all_dx) == self.getMovieLength())
        for i in range(self.getMovieLength()):
            try:
                self.setDriftCorrection(i,
                                        dx = all_dx[i],
                                        dy = all_dy[i],
                                        dz = all_dz[i])
            except saH5Py.SAH5PyException:
                pass
        
    def setDriftCorrectionXY(self, dx, dy):
        self.dx = dx
        self.dy = dy

    def setDriftCorrectionZ(self, dz):
        self.dz = dz
        
    def setFrameRange(self, fmin, fmax):
        self.fmin = fmin
        self.fmax = fmax


class SAH5DriftCorrectionTest(SAH5DriftCorrection):
    """
    A sub-class of SAH5PyDriftCorrection for testing purposes.
    """
    def locsInFrameRangeIterator(self, start, stop, fields):
        for i in range(start, stop):
            locs = self.getLocalizationsInFrame(i,
                                                drift_corrected = True,
                                                fields = fields)
            yield locs


def interpolateData(xvals, yvals, film_l):
    """
    Interpolate drift data to the length of the film.
    """

    final_drift = numpy.zeros(film_l)
    
    # Use polyfit for extrapolation at the end points.
    pe = numpy.poly1d(numpy.polyfit(xvals[0:2], yvals[0:2], 1))
    for i in range(int(xvals[0])):
        final_drift[i] = pe(i)

    pe = numpy.poly1d(numpy.polyfit(xvals[-2:], yvals[-2:], 1))
    for i in range(int(xvals[-1]), film_l):
        final_drift[i] = pe(i)        

    # Create linear spline for interpolation.
    sp = scipy.interpolate.interp1d(xvals, yvals, kind = "linear")

    # Interpolate.
    i = int(xvals[0])
    while (i <= int(xvals[-1])):
        final_drift[i] = sp(i)
        i += 1

    return final_drift

def saveDriftData(filename, fdx, fdy, fdz):
    """
    Save the x,y and z drift data to a file.
    """
    frames = numpy.arange(fdx.size) + 1
    numpy.savetxt(filename,
                  numpy.column_stack((frames,
                                      -fdx, 
                                      -fdy, 
                                      -fdz)),
                  fmt = "%d\t%.3f\t%.3f\t%.3f")

    
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
