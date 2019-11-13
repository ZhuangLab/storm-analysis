#!/usr/bin/env python
"""
Pure Python code for doing ADMM in 3D.

Hazen 11/19
"""
import numpy

import storm_analysis.sa_library.recenter_psf as recenterPSF

import storm_analysis.admm.admm_math as admmMath


class ADMM(object):

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
        self.image = None
        self.rho = rho
        self.shape = psfs.shape
        self.u = None
        self.x = None
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

    def cleanup(self):
        pass

    def getXVector(self):
        return self.x
    
    def iterate(self, l_term):

        # X update.
        self.x = admmMath.multiplyMatVec(self.G_inv, self.Atb + self.rho * (self.z - self.u))

        # Z update.
        self.z = self.x + self.u - l_term/self.rho
        self.z[(self.z < 0.0)] = 0.0

        # U update.
        self.u = self.u + self.x - self.z

    def l1Error(self):
        return numpy.sum(numpy.abs(self.x))

    def l2Error(self):
        Ax = admmMath.multiplyMatVec(self.A, self.x)
        assert(Ax.shape[2] == 1)
        return numpy.linalg.norm(Ax[:,:,0] - self.image)

    def newImage(self, image):
        self.image = image
        self.Atb = admmMath.multiplyMatVec(self.At, image)

        self.x = numpy.zeros(self.shape)
        self.z = numpy.zeros(self.shape)
        self.u = numpy.zeros(self.shape)

    
