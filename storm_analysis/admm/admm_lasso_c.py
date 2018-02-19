#!/usr/bin/env python
"""
Python interface to ADMM Lasso C library.

Hazen 02/18
"""
import ctypes
import numpy
from numpy.ctypeslib import ndpointer

import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.recenter_psf as recenterPSF


admm_lasso = loadclib.loadCLibrary("storm_analysis.admm", "admm_lasso")

# C interface definition.
admm_lasso.getXVector.argtypes = [ctypes.c_void_p, ndpointer(dtype=numpy.float64)]

admm_lasso.initialize.argtypes = [ndpointer(dtype=numpy.float64),
                                  ctypes.c_double,
                                  ctypes.c_int,
                                  ctypes.c_int]
admm_lasso.initialize.restype = ctypes.c_void_p

admm_lasso.iterate.argtypes = [ctypes.c_void_p,
                               ctypes.c_double,
                               ctypes.c_int]

admm_lasso.newImage.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype=numpy.float64)]


class ADMMLassoException(Exception):
    pass


class ADMMLasso(object):
    """
    ADMM Solver class.
    """
    def __init__(self, psfs, rho):
        super(ADMMLasso, self).__init__()

        self.shape = psfs.shape

        c_psf = numpy.ascontiguousarray(recenterPSF.recenterPSF(psfs[:,:,0]), dtype = numpy.float64)
        self.c_admm_lasso = admm_lasso.initialize(c_psf, rho, self.shape[0], self.shape[1])

    def cleanup(self):
        admm_lasso.cleanup(self.c_admm_lasso)
        self.c_admm_lasso = None

    def getXVector(self):
        c_xvec = numpy.ascontiguousarray(numpy.zeros(self.shape), dtype=numpy.float64)
        admm_lasso.getXVector(self.c_admm_lasso, c_xvec)
        return c_xvec

    def iterate(self, a_lambda, pos_only = True):
        admm_lasso.iterate(self.c_admm_lasso, a_lambda, pos_only)

    def l2Error(self):
        #
        # FIXME: Add L2 error calculation.
        #
        return "NA"
        
    def newImage(self, image):
        if (image.shape[0] != self.shape[0]) or (image.shape[1] != self.shape[1]):
            raise ADMMLassoException("image shape does not match psf shape " + " ".join([xsize, self.xsize, ysize, self.ysize]))

        c_image = numpy.ascontiguousarray(image, dtype=numpy.float64)
        admm_lasso.newImage(self.c_admm_lasso, c_image)


#
# The MIT License
#
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
