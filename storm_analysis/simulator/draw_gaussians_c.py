#!/usr/bin/env python
"""
Draws guassians onto a user supplied image.

Hazen 01/16
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.loadclib as loadclib

drawgauss = loadclib.loadCLibrary("storm_analysis.simulator", "draw_gaussians")

drawgauss.drawGaussians.argtypes = [ndpointer(dtype = numpy.float64),
                                    ndpointer(dtype = numpy.float64),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]

def cDrawGaussians(image, objects, resolution):
    c_image = numpy.ascontiguousarray(image, dtype = numpy.float64)
    c_objects = numpy.ascontiguousarray(objects, dtype = numpy.float64)
    drawgauss.drawGaussians(c_image,
                            c_objects,
                            c_image.shape[1],
                            c_image.shape[0],
                            objects.shape[0],
                            resolution)
    return c_image

def drawGaussians(size, objects, background = 0.0, res = 5):
    image = background * numpy.ones((size[0],size[1]))
    image = cDrawGaussians(image, objects, res)
    return image

def drawGaussiansXY(size, x, y, height = 1.0, sigma = 1.0, background = 0.0):
    image = background * numpy.ones((size[0],size[1]))
    np = x.shape[0]
    h = numpy.ones(np) * float(height)
    sx = numpy.ones(np) * float(sigma)
    sy = numpy.ones(np) * float(sigma)
    objects = numpy.concatenate((x[:,None],
                                 y[:,None],
                                 h[:,None],
                                 sx[:,None],
                                 sy[:,None]),
                                axis = 1)
    image = cDrawGaussians(image, objects, 5)
    return image

def drawGaussiansXYOnImage(image, x, y, height = 1.0, sigma = 1.0):
    np = x.shape[0]
    h = numpy.ones(np) * float(height)
    sx = numpy.ones(np) * float(sigma)
    sy = numpy.ones(np) * float(sigma)
    objects = numpy.concatenate((x[:,None],
                                 y[:,None],
                                 h[:,None],
                                 sx[:,None],
                                 sy[:,None]),
                                axis = 1)
    image = cDrawGaussians(image, objects, 1)
    return image

if (__name__ == "__main__"):
    import sa_library.daxwriter as daxwriter

    x = 4.0*numpy.arange(10) + 10
    y = 4.0*numpy.arange(10) + 5

    image = drawGaussiansXY([100,100], x, y, height = 100.0)
    daxwriter.singleFrameDax("dg_test.dax", image)


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
