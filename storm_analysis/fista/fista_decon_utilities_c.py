#!/usr/bin/python
#
# Simple Python interface to fista_decon_utilities.c
#
# Hazen 1/16
#

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import sa_library.loadclib as loadclib

fd_util = loadclib.loadCLibrary(os.path.dirname(__file__), "fista_decon_utilities")

# C interface definition
fd_util.label.argtypes = [ndpointer(dtype=numpy.float64),
                          ndpointer(dtype=numpy.int32),
                          ctypes.c_double,
                          ctypes.c_int,
                          ctypes.c_int,
                          ctypes.c_int]
fd_util.label.restype = ctypes.c_int
fd_util.moments.argtypes = [ndpointer(dtype=numpy.float64),
                            ndpointer(dtype=numpy.float64),
                            ndpointer(dtype=numpy.int32),
                            ctypes.c_int,
                            ctypes.c_int,
                            ctypes.c_int,
                            ctypes.c_int]


# This is just a convenience wrapper.
def getPeaks(image, threshold, margin):
    [labeled, counts] = label(image, threshold, margin)
    return moments(image, labeled, counts)


# Return a labeled version of a FISTA deconvolved 3D image.
def label(image, threshold, margin):
    image_c = numpy.ascontiguousarray(image, dtype = numpy.float64)
    labels_c = numpy.ascontiguousarray(numpy.zeros(image.shape, dtype = numpy.int32))
    counts = fd_util.label(image_c,
                           labels_c,
                           threshold,
                           margin,
                           image_c.shape[0],
                           image_c.shape[1],
                           image_c.shape[2])
    return [labels_c, counts]


# Given image and labels, return the sum and center of mass
# of each labeled part of the image.
def moments(image, labels, counts):
    image_c = numpy.ascontiguousarray(image, dtype = numpy.float64)
    labels_c = numpy.ascontiguousarray(labels, dtype = numpy.int32)
    peaks_c = numpy.ascontiguousarray(numpy.zeros((counts, 4), dtype = numpy.float64))
    fd_util.moments(image_c,
                    peaks_c,
                    labels_c,
                    counts,
                    image_c.shape[0],
                    image_c.shape[1],
                    image_c.shape[2])
    return peaks_c


if (__name__ == "__main__"):
    
    import sa_library.daxwriter as daxwriter

    # Create test image with some non-zero data.
    test_image = numpy.zeros((10,10,3))

    test_image[0,0,0] = 1.0
    
    test_image[5,5,0] = 1.0
    test_image[5,5,1] = 2.0

    test_image[7,6,1] = 1.0
    test_image[7,7,1] = 2.0
    test_image[7,8,1] = 1.0

    [labels, counts] = label(test_image, 0.1, 0)
    #print "1"
    #peaks = moments(test_image, labels, counts)
    #print "2"

    peaks = getPeaks(test_image, 0.1, 0)

    print peaks
    
    labels_image = daxwriter.DaxWriter("fd_util_test.dax", test_image.shape[0], test_image.shape[1])
    for i in range(labels.shape[2]):
        labels_image.addFrame(labels[:,:,i])
    labels_image.close()

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
