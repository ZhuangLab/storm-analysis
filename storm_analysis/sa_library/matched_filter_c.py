#!/usr/bin/env python
#
# Simply Python interface to matched_filter.c
#
# Hazen 3/16
#

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.recenter_psf as recenterPSF

m_filter = loadclib.loadCLibrary("storm_analysis.sa_library", "_matched_filter")

m_filter.cleanup.argtypes = [ctypes.c_void_p]
m_filter.convolve.argtypes = [ctypes.c_void_p,
                              ndpointer(dtype = numpy.float64),
                              ndpointer(dtype = numpy.float64)]
m_filter.initialize.argtypes = [ndpointer(dtype = numpy.float64),
                                ctypes.c_int,
                                ctypes.c_int,
                                ctypes.c_int]
m_filter.initialize.restype = ctypes.c_void_p


class MatchedFilterException(Exception):

    def __init__(self, message):
        Exception.__init__(self, message)

        
class MatchedFilter(object):

    def __init__(self, psf, estimate_fft_plan = False):
        """
        If you are only going to use this object on a few images using 
        estimate_fft_plan = True is a good idea as the initialization
        will go a lot faster, particularly for large images.
        """
        self.psf_shape = psf.shape

        rc_psf = recenterPSF.recenterPSF(psf)

        self.mfilter = m_filter.initialize(rc_psf, rc_psf.shape[0], rc_psf.shape[1], int(estimate_fft_plan))

    def cleanup(self):
        m_filter.cleanup(self.mfilter)

    def convolve(self, image):
        if (image.shape[0] != self.psf_shape[0]) or (image.shape[1] != self.psf_shape[1]):
            raise MatchedFilterException("Image shape must match psf shape! " + str(image.shape) + " != " + str(self.psf_shape))
        
        image = numpy.ascontiguousarray(image, dtype = numpy.float64)
        result = numpy.zeros(self.psf_shape, dtype = numpy.float64)
        m_filter.convolve(self.mfilter, image, result)

        return result
        

if (__name__ == "__main__"):
    
    import tifffile

    import simulator.drawgaussians as dg

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
        tif.save(result1.astype(numpy.uint16))
        tif.save(result2.astype(numpy.uint16))


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
