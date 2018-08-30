#!/usr/bin/env python
"""
Python interface for affine_transform.c

Hazen 5/17
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

import storm_analysis.sa_library.loadclib as loadclib

a_trans = loadclib.loadCLibrary("affine_transform")

a_trans.cleanup.argtypes = [ctypes.c_void_p]
a_trans.transform.argtypes = [ctypes.c_void_p,
                              ndpointer(dtype = numpy.float64),
                              ndpointer(dtype = numpy.float64),
                              ctypes.c_int,
                              ctypes.c_int]
a_trans.initialize.argtypes = [ndpointer(dtype = numpy.float64),
                               ndpointer(dtype = numpy.float64)]
a_trans.initialize.restype = ctypes.c_void_p


class AffineTransform(object):

    def __init__(self, xt = None, yt = None, **kwds):
        """
        xt and yt are the array specifying how to transform the x and y coordinates.

        xf = xt[0] + xt[1] * xi + xt[2] * yi
        yf = yt[0] + yt[1] * xi + yt[2] * yi

        The convention here is the same as that of sa_library/multi_fit, 'x' is the
        fast axis and 'y' is the slow axis.
        """
        super().__init__(**kwds)
        self.atrans = a_trans.initialize(numpy.ascontiguousarray(xt, dtype = numpy.float64),
                                         numpy.ascontiguousarray(yt, dtype = numpy.float64))

    def cleanup(self):
        a_trans.cleanup(self.atrans)

    def transform(self, image):
        iy = image.shape[0]
        ix = image.shape[1]

        trans_image = numpy.ascontiguousarray(numpy.zeros(image.shape, dtype = numpy.float64))
        a_trans.transform(self.atrans,
                          numpy.ascontiguousarray(image, dtype = numpy.float64),
                          trans_image,
                          iy,
                          ix)
        
        return trans_image


if (__name__ == "__main__"):

    import tifffile

    # Make a grid image.
    image = numpy.zeros((300,200))
    for i in range(10):
        image[i*30+15,:] = 1.0
        image[:,i*20+10] = 1.0

    at = AffineTransform(xt = [-20.0, 1.1, 0.1],
                         yt = [-10.0, 0.1, 1.1])

    tr_image = at.transform(image)

    with tifffile.TiffWriter("transform.tif") as tf:
        tf.save(image.astype(numpy.float32))
        tf.save(tr_image.astype(numpy.float32))

    at.cleanup()
   

#
# The MIT License
#
# Copyright (c) 2017 Babcock Lab, Harvard University
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
