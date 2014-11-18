#!/usr/bin/python
#
# Simple Python interface to the sCMOS image processing C library.
# This is based on the algorithms described in:
#
# "Video-rate nanoscopy using sCMOS camera-specific single-molecule localization algorithms"
# F. Huang et al. Nature Methods, 10, p653-658.
#
# Hazen 10/13
#

import ctypes
import math
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import sa_library.loadclib as loadclib

slib = loadclib.loadCLibrary(os.path.dirname(__file__), "scmos_utilities")

# C interface definition.
slib.deregularize.argtypes = [ndpointer(dtype=numpy.float64),
                              ndpointer(dtype=numpy.float64),
                              ndpointer(dtype=numpy.float64),
                              ndpointer(dtype=numpy.float64),
                              ndpointer(dtype=numpy.float64),
                              ctypes.c_int]

slib.regularize.argtypes = [ndpointer(dtype=numpy.float64),
                            ndpointer(dtype=numpy.float64),
                            ndpointer(dtype=numpy.float64),
                            ndpointer(dtype=numpy.float64),
                            ndpointer(dtype=numpy.float64),
                            ctypes.c_int]

slib.smooth.argtypes = [ndpointer(dtype=numpy.float64),
                        ndpointer(dtype=numpy.float64),
                        ndpointer(dtype=numpy.float64),
                        ndpointer(dtype=numpy.float64),
                        ndpointer(dtype=numpy.float64),
                        ctypes.c_int,
                        ctypes.c_int,
                        ctypes.c_int,
                        ctypes.c_int]

#
# Image regularizer class.
#
# See Section 3.2 in the supplement, "Single-Particle Localization using MLEsCMOS.
#
class Regularizer:
    def __init__(self, camera_offset, camera_variance, camera_gain):
        self.c_offset = numpy.ascontiguousarray(camera_offset, dtype = numpy.float64)
        self.c_gain = numpy.ascontiguousarray(camera_gain, dtype = numpy.float64)
        self.c_variance = numpy.ascontiguousarray(camera_variance, dtype = numpy.float64)
        self.shape = self.c_offset.shape

    def deregularizeImage(self, image):
        assert (image.size == self.c_offset.size), "Images must be the same size!"
        c_deregularized = numpy.ascontiguousarray(numpy.zeros(self.shape), dtype = numpy.float64)
        c_image = numpy.ascontiguousarray(image, dtype = numpy.float64)
        slib.deregularize(c_deregularized,
                          c_image,
                          self.c_offset,
                          self.c_variance,
                          self.c_gain,
                          image.size)
        return c_deregularized

    def regularizeImage(self, image):
        assert (image.size == self.c_offset.size), "Images must be the same size!"
        c_regularized = numpy.ascontiguousarray(numpy.zeros(self.shape), dtype = numpy.float64)
        c_image = numpy.ascontiguousarray(image, dtype = numpy.float64)
        slib.regularize(c_regularized,
                        c_image,
                        self.c_offset,
                        self.c_variance,
                        self.c_gain,
                        image.size)
        return c_regularized

    def normalizeImage(self, image):
        return image/self.c_gain - self.c_offset

#
# Image smoother class.
#
# See Section 3.1 in the supplement, "Image Segmentation".
#
class Smoother:

    def __init__(self, camera_offset, camera_variance, camera_gain):
        self.c_offset = numpy.ascontiguousarray(camera_offset, dtype = numpy.float64)
        self.c_inv_gain = numpy.ascontiguousarray(1.0/camera_gain, dtype = numpy.float64)
        self.c_inv_variance = numpy.ascontiguousarray(1.0/camera_variance, dtype = numpy.float64)

        self.shape = camera_offset.shape
        self.x_size = self.shape[0]
        self.y_size = self.shape[1]

    def smoothImage(self, image, sigma_psf):
        assert (image.size == self.c_offset.size), "Images must be the same size!"
        c_smoothed = numpy.ascontiguousarray(numpy.zeros(self.shape), dtype = numpy.float64)
        c_image = numpy.ascontiguousarray(image, dtype = numpy.float64)
        c_sigma_psf = int(round(sigma_psf))
        slib.smooth(c_smoothed,
                    c_image,
                    self.c_offset,
                    self.c_inv_variance,
                    self.c_inv_gain,
                    self.x_size,
                    self.y_size,
                    c_sigma_psf,
                    2*c_sigma_psf)
        return c_smoothed


# Testing.
if __name__ == "__main__":

    xsize = 3
    ysize = 3

    gain = numpy.ones((xsize,ysize))
    offset = 10.0*numpy.ones((xsize,ysize))
    variance = 0.5*numpy.ones((xsize,ysize))

    regf = Regularizer(offset, variance, gain)

    image = numpy.ones((xsize,ysize))
    print image

    r_image = regf.regularizeImage(image)
    print r_image

    rr_image = regf.deregularizeImage(r_image)
    print rr_image

#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
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
