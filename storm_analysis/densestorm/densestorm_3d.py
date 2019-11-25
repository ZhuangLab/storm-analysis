#!/usr/bin/env python
"""
Pure Python code for Poisson sparse deconvolution following 
3DenseSTORM.

As described in:
Ovesny et al., "High density 3D localization microscopy using
sparse support recovery", Optics Express, 2014.

Hazen 11/19
"""
import numpy

import storm_analysis.sa_library.cs_algorithm as csAlgorithm
import storm_analysis.sa_library.recenter_psf as recenterPSF

import storm_analysis.admm.admm_math as admmMath


class DenseSTORM(csAlgorithm.CSAlgorithm):

    def __init__(self, psfs, beta, eta, micro, **kwds):
        """
        psfs is an array of psfs for different image planes (nx, ny, nz).
        They must be the same size as the image that will get analyzed.
        """        
        super(DenseSTORM, self).__init__(**kwds)

        #
        # Note: Use A instead of H, but the notation should
        #       otherwise be consistent with Ovesny et al.
        #
        # _h = hat,_t = tilde
        #
        self.A = None
        self.At = None
        self.b_vec = None
        self.beta = beta
        self.d_vec = None
        self.e_vec = None
        self.eta = eta
        self.G_inv = None
        self.micro = micro
        self.shape = psfs.shape
        self.w_vec = None
        self.x_vec_h = None
        self.x_vec_t = None
        self.y_vec = None
        self.y_vec_h = None
        
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
            G[i,i] += admmMath.identityMatrix(G.getMatrixShape(), scale = self.micro)

        [L,D,U] = admmMath.lduG(G)

        L_inv = admmMath.invL(L)
        D_inv = admmMath.invD(D)
        U_inv = admmMath.invU(U)
        
        self.G_inv = admmMath.multiplyMatMat(U_inv, admmMath.multiplyMatMat(D_inv, L_inv))
    
    def getAx(self):
        Ax = admmMath.multiplyMatVec(self.A, self.x_vec_h)
        assert(Ax.shape[2] == 1)
        return Ax[:,:,0]

    def getXVector(self):
        return self.x_vec_h
    
    def iterate(self):

        # Equation 6.
        At_ydb = admmMath.multiplyMatVec(self.At, self.y_vec_h + self.d_vec - self.b_vec)
        self.x_vec_h = admmMath.multiplyMatVec(self.G_inv, At_ydb + self.micro * (self.x_vec_t + self.e_vec))

        # Equation 7.
        Ax = self.getAx()
        m_Ax_bd = self.eta*(Ax + self.b_vec - self.d_vec)
        t1 = 1.0 - m_Ax_bd
        self.y_vec_h = (-t1 + numpy.sqrt(t1*t1 + 4.0*self.eta*self.y_vec))/(2.0*self.eta)

        # Equation 8.
        self.x_vec_t = self.x_vec_h - self.e_vec - self.w_vec
        self.x_vec_t[(self.x_vec_t < 0.0)] = 0.0

        # Equation 9.
        self.d_vec = self.d_vec - (Ax + self.b_vec - self.y_vec_h)

        # Equation 10.
        self.e_vec = self.e_vec - (self.x_vec_h - self.x_vec_t)

    def newImage(self, image, background):
        """
        image - The image to deconvolve, includes the background estimate.
        background - An estimate of the background in the image.
        """
        # These are used by the super-class.
        self.image = image - background
        self.weights = 1.0/image
        
        # Fixed vectors.
        self.y_vec = image
        self.b_vec = background

        t1 = self.beta*numpy.sqrt(background)/self.micro
        self.w_vec = numpy.zeros(self.shape)
        for i in range(self.w_vec.shape[2]):
            self.w_vec[:,:,i] = t1
        
        # Updated vectors.
        self.e_vec = numpy.zeros(self.shape)
        self.x_vec_h = numpy.zeros(self.shape)
        self.x_vec_t = numpy.zeros(self.shape)

        self.d_vec = numpy.zeros((self.shape[0], self.shape[1]))
        self.y_vec_h = numpy.copy(self.y_vec)
