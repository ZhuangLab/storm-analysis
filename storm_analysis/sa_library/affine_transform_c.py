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

a_trans = loadclib.loadCLibrary("storm_analysis.sa_library", "affine_transform")

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
        """
        super().__init__(**kwds)

        self.atrans = a_trans.initialize(numpy.ascontiguousarray(xt, dtype = numpy.float64),
                                         numpy.ascontiguousarray(yt, dtype = numpy.float64))

    def cleanup(self):
        a_trans.cleanup(self.atrans)

    def transform(self, image):
        ix = image.size[0]
        iy = image.size[1]

        trans_image = numpy.ascontiguousarray(numpy.zeros(image.shape, dtype = numpy.float64))
        a_trans.transform(self.atrans,
                          numpy.ascontiguousarray(image, dtype = numpy.float64),
                          trans_image,
                          ix,
                          iy)
        
        return trans_image
