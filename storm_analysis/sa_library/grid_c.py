#!/usr/bin/env python
"""
Simple Python interface to grid.c. This is a somewhat faster
but less flexible approach to creating 2D and 3D histograms
than using the built-in numpy function numpy.histogramdd().

Hazen 12/11
"""

from ctypes import *
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.loadclib as loadclib

grid = loadclib.loadCLibrary("grid")

# Function specifications
grid.grid2D.argtypes = [ndpointer(dtype=numpy.int32),
                        ndpointer(dtype=numpy.int32),
                        ndpointer(dtype=numpy.int32),
                        c_int,
                        c_int,
                        c_int]

grid.grid3D.argtypes = [ndpointer(dtype=numpy.int32),
                        ndpointer(dtype=numpy.int32),
                        ndpointer(dtype=numpy.int32),
                        ndpointer(dtype=numpy.int32),
                        c_int,
                        c_int,
                        c_int,
                        c_int]

def grid2D(x,y,image):
    """
    Grid in 2D.
    """
    assert (image.dtype == numpy.int32)
    assert (image.flags['C_CONTIGUOUS']), "Image is not C contiguous."
    assert (len(image.shape) == 2)
    
    c_x = numpy.ascontiguousarray(x).astype(numpy.int32)
    c_y = numpy.ascontiguousarray(y).astype(numpy.int32)
    grid.grid2D(image,
                c_x,
                c_y,
                image.shape[0],
                image.shape[1],
                c_x.size)

def grid3D(x,y,z,image):
    """
    Grid in 3D.
    """
    assert (image.dtype == numpy.int32)
    assert (image.flags['C_CONTIGUOUS']), "Image is not C contiguous."
    assert (len(image.shape) == 3)
    
    c_x = numpy.ascontiguousarray(x).astype(numpy.int32)
    c_y = numpy.ascontiguousarray(y).astype(numpy.int32)
    c_z = numpy.ascontiguousarray(z).astype(numpy.int32)
    grid.grid3D(image,
                c_x,
                c_y,
                c_z,
                image.shape[0],
                image.shape[1],
                image.shape[2],
                c_x.size)


#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
