#!/usr/bin/env python
"""
Python interface to cs_decon_utilities.c

Hazen 2/18
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.loadclib as loadclib

cs_util = loadclib.loadCLibrary("cs_decon_utilities")

# C interface definition
cs_util.label.argtypes = [ndpointer(dtype=numpy.float64),
                          ndpointer(dtype=numpy.int32),
                          ctypes.c_double,
                          ctypes.c_int,
                          ctypes.c_int,
                          ctypes.c_int]
cs_util.label.restype = ctypes.c_int
cs_util.moments.argtypes = [ndpointer(dtype=numpy.float64),
                            ndpointer(dtype=numpy.float64),
                            ndpointer(dtype=numpy.int32),
                            ctypes.c_int,
                            ctypes.c_int,
                            ctypes.c_int,
                            ctypes.c_int]


def getPeaks(image, threshold, margin):
    """
    This is just a convenience wrapper.
    """
    [labeled, counts] = label(image, threshold, margin)
    return moments(image, labeled, counts)


def label(image, threshold, margin):
    """
    Return a labeled version of a CS deconvolved 3D image.
    """
    #
    # Note: Add 1 to margin to match what we do in ia_utilities_c.py.
    #
    image_c = numpy.ascontiguousarray(image, dtype = numpy.float64)
    labels_c = numpy.ascontiguousarray(numpy.zeros(image.shape, dtype = numpy.int32))
    counts = cs_util.label(image_c,
                           labels_c,
                           threshold,
                           margin+1, 
                           image_c.shape[0],
                           image_c.shape[1],
                           image_c.shape[2])
    return [labels_c, counts]


def moments(image, labels, counts):
    """
    Given image and labels, return the sum and center of mass
    of each labeled part of the image.
    """
    image_c = numpy.ascontiguousarray(image, dtype = numpy.float64)
    labels_c = numpy.ascontiguousarray(labels, dtype = numpy.int32)
    peaks_c = numpy.ascontiguousarray(numpy.zeros((counts, 4), dtype = numpy.float64))
    cs_util.moments(image_c,
                    peaks_c,
                    labels_c,
                    counts,
                    image_c.shape[0],
                    image_c.shape[1],
                    image_c.shape[2])
    return peaks_c


#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
# Copyright (c) 2018 Babcock Lab, Harvard University
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
