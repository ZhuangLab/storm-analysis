#!/usr/bin/env python
"""
Python interface to ADMM Lasso C library.

Hazen 11/19
"""
import ctypes
import numpy
from numpy.ctypeslib import ndpointer

import storm_analysis.sa_library.cs_algorithm as csAlgorithm
import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.recenter_psf as recenterPSF

import storm_analysis.admm.admm_math as admmMath


admm_lasso = loadclib.loadCLibrary("admm_lasso")

# C interface definition.
admm_lasso.getAx.argtypes = [ctypes.c_void_p,
                             ndpointer(dtype=numpy.float64)]

admm_lasso.getXVector.argtypes = [ctypes.c_void_p,
                                  ndpointer(dtype=numpy.float64)]

admm_lasso.initialize2D.argtypes = [ctypes.c_double,
                                    ctypes.c_int,
                                    ctypes.c_int]
admm_lasso.initialize2D.restype = ctypes.c_void_p

admm_lasso.initialize3D.argtypes = [ctypes.c_double,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]
admm_lasso.initialize3D.restype = ctypes.c_void_p

admm_lasso.initializeA.argtypes = [ctypes.c_void_p,
                                   ndpointer(dtype=numpy.complex128),
                                   ctypes.c_int]

admm_lasso.initializeGInv.argtypes = [ctypes.c_void_p,
                                      ndpointer(dtype=numpy.complex128),
                                      ctypes.c_int]

admm_lasso.iterate.argtypes = [ctypes.c_void_p,
                               ctypes.c_double,
                               ctypes.c_int]

admm_lasso.l1Error.argtypes = [ctypes.c_void_p]
admm_lasso.l1Error.restype = ctypes.c_double

admm_lasso.l2Error.argtypes = [ctypes.c_void_p]
admm_lasso.l2Error.restype = ctypes.c_double

admm_lasso.newImage.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype=numpy.float64)]

admm_lasso.run.argtypes = [ctypes.c_void_p,
                           ctypes.c_double,
                           ctypes.c_int]


class ADMMLassoException(Exception):
    pass


class ADMMLasso(csAlgorithm.CSAlgorithm):
    """
    ADMM Solver class.
    """
    def __init__(self, psfs, rho, **kwds):
        super(ADMMLasso, self).__init__(**kwds)

        self.shape = psfs.shape

        # Initialize C library.
        self.c_admm_lasso = admm_lasso.initialize3D(rho, self.shape[0], self.shape[1], self.shape[2])

        #
        # Do the ADMM math on the Python side.
        #
        # Calculate A matrix.
        nz = self.shape[2]
        A = admmMath.Cells(nz, 1)
        for i in range(nz):
            tmp = recenterPSF.recenterPSF(psfs[:,:,i])
            A[i,0] = numpy.fft.fft2(tmp)

        # Calculate A transpose.
        At = admmMath.transpose(A)

        # Calculate AtA + rhoI inverse.
        G = admmMath.multiplyMatMat(At, A)

        for i in range(nz):
            G[i,i] += admmMath.identityMatrix(G.getMatrixShape(), scale = rho)

        [L,D,U] = admmMath.lduG(G)

        L_inv = admmMath.invL(L)
        D_inv = admmMath.invD(D)
        U_inv = admmMath.invU(U)
        
        G_inv = admmMath.multiplyMatMat(U_inv, admmMath.multiplyMatMat(D_inv, L_inv))
        
        # Initialize A and G_inv matrices in the C library.
        fft_size = int(self.shape[1]/2+1)
        for i in range(nz):

            # Remove redundant frequencies that FFTW doesn't use.
            c_A = A[i,0][:,:fft_size]
            c_A = numpy.ascontiguousarray(c_A, dtype = numpy.complex128)
            admm_lasso.initializeA(self.c_admm_lasso, c_A, i)

        for i in range(nz):
            for j in range(nz):
                
                # Remove redundant frequencies that FFTW doesn't use.
                #
                # We index (j,i) here because this is what gives us results that
                # match admm_3d (the pure Python version of 3D ADMM.
                #
                c_G_inv = G_inv[j,i][:,:fft_size]
                c_G_inv = numpy.ascontiguousarray(c_G_inv, dtype = numpy.complex128)
                admm_lasso.initializeGInv(self.c_admm_lasso, c_G_inv, i*nz + j)

    def cleanup(self):
        admm_lasso.cleanup(self.c_admm_lasso)
        self.c_admm_lasso = None

    def dwlsError(self):
        return 0.0
    
    def getAx(self):
        ax = numpy.ascontiguousarray(numpy.zeros((self.shape[0], self.shape[1]), dtype = float))
        admm_lasso.getAx(self.c_admm_lasso, ax)
        return ax
    
    def getXVector(self):
        c_xvec = numpy.ascontiguousarray(numpy.zeros(self.shape), dtype=numpy.float64)
        admm_lasso.getXVector(self.c_admm_lasso, c_xvec)
        return c_xvec

    def iterate(self, a_lambda, pos_only = False):
        admm_lasso.iterate(self.c_admm_lasso, a_lambda, pos_only)

    def l1Error(self):
        return admm_lasso.l1Error(self.c_admm_lasso)

    def l2Error(self):
        return admm_lasso.l2Error(self.c_admm_lasso)
    
    def newImage(self, image, background):
        image_no_bg = image - background
        
        if (image_no_bg.shape[0] != self.shape[0]) or (image_no_bg.shape[1] != self.shape[1]):
            raise ADMMLassoException("Image shape does not match psf shape " + " ".join([str(image_no_bg.shape), str(self.shape)]))

        c_image = numpy.ascontiguousarray(image_no_bg, dtype=numpy.float64)
        admm_lasso.newImage(self.c_admm_lasso, c_image)

    def run(self, a_lamba, iterations):
        admm_lasso.run(self.c_admm_lasso, a_lambda, iterations)


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
