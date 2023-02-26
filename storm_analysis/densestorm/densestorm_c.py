#!/usr/bin/env python
"""
Python interface to 3denseSTORM C library.

Hazen 11/19
"""
import ctypes
import numpy
from numpy.ctypeslib import ndpointer

import storm_analysis.sa_library.cs_algorithm as csAlgorithm
import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.recenter_psf as recenterPSF

import storm_analysis.admm.admm_math as admmMath


densestorm = loadclib.loadCLibrary("densestorm")

# C interface definition.
densestorm.getAx.argtypes = [ctypes.c_void_p,
                             ndpointer(dtype=numpy.float64)]

densestorm.getXVector.argtypes = [ctypes.c_void_p,
                                  ndpointer(dtype=numpy.float64),
                                  ctypes.c_int]

densestorm.initialize2D.argtypes = [ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_int,
                                    ctypes.c_int]
densestorm.initialize2D.restype = ctypes.c_void_p

densestorm.initialize3D.argtypes = [ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]
densestorm.initialize3D.restype = ctypes.c_void_p

densestorm.initializeA.argtypes = [ctypes.c_void_p,
                                   ndpointer(dtype=numpy.complex128),
                                   ctypes.c_int]

densestorm.initializeGInv.argtypes = [ctypes.c_void_p,
                                      ndpointer(dtype=numpy.complex128),
                                      ctypes.c_int]

densestorm.iterate.argtypes = [ctypes.c_void_p]

densestorm.l1Error.argtypes = [ctypes.c_void_p]
densestorm.l1Error.restype = ctypes.c_double

densestorm.l2Error.argtypes = [ctypes.c_void_p]
densestorm.l2Error.restype = ctypes.c_double

densestorm.newImage.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype=numpy.float64),
                                ndpointer(dtype=numpy.float64)]

densestorm.run.argtypes = [ctypes.c_void_p,
                           ctypes.c_int]


class DenseSTORMException(Exception):
    pass


class DenseSTORM(csAlgorithm.CSAlgorithm):
    """
    3denseSTORM Solver class.
    """
    def __init__(self, psfs, beta, eta, micro, **kwds):
        super(DenseSTORM, self).__init__(**kwds)

        self.shape = psfs.shape

        # Initialize C library.
        self.c_densestorm = densestorm.initialize3D(beta, eta, micro, self.shape[0], self.shape[1], self.shape[2])

        #
        # Do the A, G matrix math in Python.
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
            G[i,i] += admmMath.identityMatrix(G.getMatrixShape(), scale = micro)

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
            densestorm.initializeA(self.c_densestorm, c_A, i)

        for i in range(nz):
            for j in range(nz):
                
                # Remove redundant frequencies that FFTW doesn't use.
                #
                # We index (j,i) here because this is what gives us results that
                # match admm_3d (the pure Python version of 3D ADMM.
                #
                c_G_inv = G_inv[j,i][:,:fft_size]
                c_G_inv = numpy.ascontiguousarray(c_G_inv, dtype = numpy.complex128)
                densestorm.initializeGInv(self.c_densestorm, c_G_inv, i*nz + j)

    def cleanup(self):
        densestorm.cleanup(self.c_densestorm)
        self.c_densestorm = None

    def dwlsError(self):
        return 0.0
    
    def getAx(self):
        ax = numpy.ascontiguousarray(numpy.zeros((self.shape[0], self.shape[1]), dtype = float))
        densestorm.getAx(self.c_densestorm, ax)
        return ax
    
    def getXVector(self, compressed = True):
        c_xvec = numpy.ascontiguousarray(numpy.zeros(self.shape), dtype=numpy.float64)
        densestorm.getXVector(self.c_densestorm, c_xvec, compressed)
        return c_xvec

    def iterate(self):
        densestorm.iterate(self.c_densestorm)

    def l1Error(self):
        return densestorm.l1Error(self.c_densestorm)

    def l2Error(self):
        return densestorm.l2Error(self.c_densestorm)
    
    def newImage(self, image, background):
        
        if (image.shape[0] != self.shape[0]) or (image.shape[1] != self.shape[1]):
            msg = "Image shape does not match psf shape " + " ".join([str(image.shape), str(self.shape)])
            raise DenseSTORMException(msg)

        if (background.shape[0] != self.shape[0]) or (background.shape[1] != self.shape[1]):
            msg = "Background shape does not match psf shape " + " ".join([str(background.shape), str(self.shape)])
            raise DenseSTORMException(msg)

        c_background = numpy.ascontiguousarray(background, dtype=numpy.float64)
        c_image = numpy.ascontiguousarray(image, dtype=numpy.float64)
        densestorm.newImage(self.c_densestorm, c_image, c_background)

    def run(self, iterations):
        densestorm.run(self.c_densestorm, iterations)


#
# The MIT License
#
# Copyright (c) 2019 Babcock Lab, Harvard University
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
