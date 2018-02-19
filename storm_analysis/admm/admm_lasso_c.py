#!/usr/bin/python
#
# Simple Python interface to ADMM Lasso C library
#
# Hazen 04/14
#

from ctypes import *
import math
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import sa_library.loadclib as loadclib

admm_lasso = loadclib.loadCLibrary(os.path.dirname(__file__), "admm_lasso")

# C interface definition.
admm_lasso.getXVector.argtypes = [ndpointer(dtype=numpy.float64)]

admm_lasso.initialize.argtypes = [ndpointer(dtype=numpy.float64),
                                  c_double,
                                  c_int]
admm_lasso.iterate.argtypes = [c_double,
                               c_int]
admm_lasso.newImage.argtypes = [ndpointer(dtype=numpy.float64)]


class ADMMLassoException(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

# Solver class.
class ADMMLasso(object):

    def __init__(self, psf, rho):
        self.initialized = True
        self.iterations = 0

        # Adjust PSF for FFT use.
        [self.xsize, self.ysize] = psf.shape
        if (self.xsize != self.ysize):
            raise ADMMLassoException("psf (and image) must be square.")

        admm_psf = numpy.roll(psf, self.xsize/2, axis = 0)
        admm_psf = numpy.roll(admm_psf, self.ysize/2, axis = 1)

        c_admm_psf = numpy.ascontiguousarray(admm_psf, dtype=numpy.float64)
        admm_lasso.initialize(c_admm_psf,
                              rho,
                              self.xsize)

    def cleanup(self):
        admm_lasso.cleanup()
        self.initialized = False

    def getXVector(self):
        xvec = numpy.zeros((self.xsize, self.ysize))
        c_xvec = numpy.ascontiguousarray(xvec, dtype=numpy.float64)
        admm_lasso.getXVector(c_xvec)
        return c_xvec

    def iterate(self, a_lambda, pos_only = True):
        c_pos_only = 0
        if pos_only:
            c_pos_only = 1
        admm_lasso.iterate(a_lambda, c_pos_only)
        self.iterations += 1
        
    def newImage(self, image):
        [xsize, ysize] = image.shape
        if (xsize != self.xsize) or (ysize != self.ysize):
            raise ADMMLassoException("image shape does not match psf shape " + " ".join([xsize, self.xsize, ysize, self.ysize]))

        c_image = numpy.ascontiguousarray(image, dtype=numpy.float64)
        admm_lasso.newImage(c_image)
        self.iterations = 0


#
# The MIT License
#
# Copyright (c) 2014 Zhuang Lab, Harvard University
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
