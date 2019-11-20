#!/usr/bin/env python
"""
Pure Python code for doing ADMM in 3D.

minimize 1/2*|| Ax - b ||_2^2 + \lambda || x ||_1

As described in:
Boyd et al., "Distributed Optimization and Statistical Learning 
via the Alternating Direction Method of Multipliers", Foundations
and Trends in Machine Learning, 2010.

Hazen 11/19
"""
import numpy

import storm_analysis.sa_library.cs_algorithm as csAlgorithm
import storm_analysis.sa_library.recenter_psf as recenterPSF

import storm_analysis.admm.admm_math as admmMath


class ADMM(csAlgorithm.CSAlgorithm):

    def __init__(self, psfs, rho, **kwds):
        """
        psfs is an array of psfs for different image planes (nx, ny, nz).
        They must be the same size as the image that will get analyzed.
        """        
        super(ADMM, self).__init__(**kwds)

        self.A = None
        self.At = None
        self.Atb = None
        self.G_inv = None
        self.rho = rho
        self.shape = psfs.shape
        self.u = None
        self.z = None

        # Calculate A matrix.
        nz = self.shape[2]
        self.A = admmMath.Cells(nz, 1)
        for i in range(nz):
            tmp = recenterPSF.recenterPSF(psfs[:,:,i])
            self.A[i,0] = numpy.fft.fft2(tmp)

        # Calculate A transpose.
        self.At = admmMath.transpose(self.A)

        # Calculate AtA + rhoI inverse.
        G = admmMath.multiplyMatMat(self.At, self.A)

        for i in range(nz):
            G[i,i] += admmMath.identityMatrix(G.getMatrixShape(), scale = rho)

        [L,D,U] = admmMath.lduG(G)

        L_inv = admmMath.invL(L)
        D_inv = admmMath.invD(D)
        U_inv = admmMath.invU(U)
        
        self.G_inv = admmMath.multiplyMatMat(U_inv, admmMath.multiplyMatMat(D_inv, L_inv))
    
    def getAx(self):
        Ax = admmMath.multiplyMatVec(self.A, self.x)
        assert(Ax.shape[2] == 1)
        return Ax[:,:,0]
    
    def iterate(self, l_term):

        # X update.
        self.x = admmMath.multiplyMatVec(self.G_inv, self.Atb + self.rho * (self.z - self.u))

        # Z update.
        self.z = self.shrink(self.x + self.u, l_term/self.rho)

        # U update.
        self.u = self.u + self.x - self.z

    def newImage(self, image, background):
        """
        image - The image to deconvolve, includes the background estimate.
        background - An estimate of the background in the image.
        """        
        self.image = image - background
        self.weights = 1.0/image
        
        self.Atb = admmMath.multiplyMatVec(self.At, self.image)

        self.x = numpy.zeros(self.shape)
        self.z = numpy.zeros(self.shape)
        self.u = numpy.zeros(self.shape)

    def shrink(self, vec, s_val):
        t1 = numpy.abs(vec) - s_val
        mask = (t1 <= 0.0)
        t1[mask] = 0.0
        
        return t1 * numpy.sign(vec)
