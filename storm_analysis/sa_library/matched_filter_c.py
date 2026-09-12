#!/usr/bin/env python
"""
Simple Python interface to matched_filter.c

Hazen 3/16
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.recenter_psf as recenterPSF

m_filter = loadclib.loadCLibrary("matched_filter")

m_filter.cleanup.argtypes = [ctypes.c_void_p]

m_filter.convolve.argtypes = [ctypes.c_void_p,
                              ndpointer(dtype = numpy.float64),
                              ndpointer(dtype = numpy.float64)]

m_filter.convolveMemo.argtypes = [ctypes.c_void_p,
                                  ndpointer(dtype = numpy.float64),
                                  ndpointer(dtype = numpy.float64)]

m_filter.initialize.argtypes = [ndpointer(dtype = numpy.float64),
                                ctypes.c_double,
                                ctypes.c_int,
                                ctypes.c_int,
                                ctypes.c_int]

m_filter.initialize.restype = ctypes.c_void_p

m_filter.cleanupBank.argtypes = [ctypes.c_void_p]

m_filter.convolveBank.argtypes = [ctypes.c_void_p,
                                  ndpointer(dtype = numpy.float64),
                                  ndpointer(dtype = numpy.float64)]

m_filter.initializeBank.argtypes = [ndpointer(dtype = numpy.float64),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]

m_filter.initializeBank.restype = ctypes.c_void_p


class MatchedFilterException(Exception):

    def __init__(self, message):
        Exception.__init__(self, message)

        
class MatchedFilter(object):

    def __init__(self, psf, fftw_estimate = False, memoize = False, max_diff = 0.1):
        """
        If you are only going to use this object on a few images using 
        fftw_estimate = True is a good idea as the initialization
        will go a lot faster, particularly for large images.

        If you think that you may being repeatedly asking it to convolve the same
        image or almost the same image then using memoization might be a good idea.
        """
        self.memoize = memoize
        
        self.psf_shape = psf.shape

        rc_psf = recenterPSF.recenterPSF(psf)

        if self.memoize:
            self.mfilter = m_filter.initialize(rc_psf,
                                               max_diff,
                                               rc_psf.shape[0],
                                               rc_psf.shape[1],
                                               int(fftw_estimate))
        else:
            self.mfilter = m_filter.initialize(rc_psf,
                                               0.0,
                                               rc_psf.shape[0],
                                               rc_psf.shape[1],
                                               int(fftw_estimate))

    def cleanup(self):
        m_filter.cleanup(self.mfilter)
        self.mfilter = None

    def convolve(self, image):
        if (image.shape[0] != self.psf_shape[0]) or (image.shape[1] != self.psf_shape[1]):
            raise MatchedFilterException("Image shape must match psf shape! " + str(image.shape) + " != " + str(self.psf_shape))
        
        image = numpy.ascontiguousarray(image, dtype = numpy.float64)
        result = numpy.zeros(self.psf_shape, dtype = numpy.float64)
        if self.memoize:
            m_filter.convolveMemo(self.mfilter, image, result)
        else:
            m_filter.convolve(self.mfilter, image, result)

        return result


class MatchedFilterBank(object):
    """
    Convolve one image with several psfs.

    This exists because applying N psfs to the same image with N MatchedFilter
    objects computes the same forward FFT of that image N times. Here it is
    computed once and reused, so the transform count per call goes from 2N to
    N+1.

    The results are the same numbers that N MatchedFilter objects produce, in
    psf order. There is no memoization, the caller is expected to be applying
    the bank to an image that has changed.
    """
    def __init__(self, psfs, fftw_estimate = False):
        """
        psfs - A list of psfs, or an (n, x, y) array. All the same shape.
        """
        psfs = numpy.asarray(psfs, dtype = numpy.float64)
        if (psfs.ndim != 3):
            raise MatchedFilterException("psfs must be (n, x, y), got " + str(psfs.shape))

        self.n_filters = psfs.shape[0]
        self.psf_shape = (psfs.shape[1], psfs.shape[2])

        rc_psfs = numpy.zeros(psfs.shape)
        for i in range(self.n_filters):
            rc_psfs[i,:,:] = recenterPSF.recenterPSF(psfs[i,:,:])
        rc_psfs = numpy.ascontiguousarray(rc_psfs, dtype = numpy.float64)

        self.mfilter = m_filter.initializeBank(rc_psfs,
                                               self.n_filters,
                                               self.psf_shape[0],
                                               self.psf_shape[1],
                                               int(fftw_estimate))

    def cleanup(self):
        m_filter.cleanupBank(self.mfilter)
        self.mfilter = None

    def convolve(self, image):
        """
        Returns an (n_filters, x, y) array, one convolution per psf.
        """
        if (image.shape[0] != self.psf_shape[0]) or (image.shape[1] != self.psf_shape[1]):
            raise MatchedFilterException("Image shape must match psf shape! " + str(image.shape) + " != " + str(self.psf_shape))

        image = numpy.ascontiguousarray(image, dtype = numpy.float64)
        results = numpy.zeros((self.n_filters, self.psf_shape[0], self.psf_shape[1]),
                              dtype = numpy.float64)
        m_filter.convolveBank(self.mfilter, image, results)

        return results


if (__name__ == "__main__"):
    
    import tifffile

    import storm_analysis.simulator.draw_gaussians_c as dg

    x_size = 200
    y_size = 300

    objects = numpy.zeros((1, 5))

    # PSF of shape 1.
    objects[0,:] = [x_size/2, y_size/2, 1.0, 1.0, 1.0]
    psf1 = dg.drawGaussians((x_size, y_size), objects)
    psf1 = psf1/numpy.sum(psf1)
    flt1 = MatchedFilter(psf1)

    # PSF of shape 2.
    objects[0,:] = [x_size/2, y_size/2, 1.0, 2.0, 0.7]
    psf2 = dg.drawGaussians((x_size, y_size), objects)
    psf2 = psf2/numpy.sum(psf2)
    flt2 = MatchedFilter(psf2)

    result1 = flt1.convolve(10000.0 * psf2)
    result2 = flt2.convolve(10000.0 * psf2)

    print("Match to 1:", numpy.max(result1))
    print("Match to 2:", numpy.max(result2))

    with tifffile.TiffWriter("result.tif") as tif:
        tif.write(result1.astype(numpy.uint16))
        tif.write(result2.astype(numpy.uint16))


#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
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
