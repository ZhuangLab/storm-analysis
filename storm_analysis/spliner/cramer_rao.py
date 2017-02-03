#!/usr/bin/env python
"""
Estimate Cramer-Rao bounds for a spline based on the formulation in
this paper (specifically equation 18):

"Three dimensional single molecule localization using a phase retrieved pupil
function", Sheng Liu, Emil B. Kromann, Wesley D. Krueger, Joerg Bewersdorf, 
and Keith A. Lidke, Optics Express 2013.

Hazen 02/17
"""

import math
import numpy
import pickle

import storm_analysis.spliner.spline_to_psf as splineToPSF


# FIXME: Also handle 2D splines.

class CramerRaoException(Exception):
    pass


class CRSplineToPSF3D(splineToPSF.SplineToPSF3D):

    def getSplineVals(self, spline_method, z_value):
        """
        Return the spline values of a given method such as:

        spline.f - The PSF
        spline.dx - The derivative in x
        ...

        This is basically a specialized version of the getPSF() method.
        """

        # Calculate PSF at requested z value.
        scaled_z = self.getScaledZ(z_value)
        
        vals_size = int((self.spline_size - 1)/2)
                
        vals = numpy.zeros((vals_size, vals_size))
        if((vals_size%2) == 0):
            for x in range(vals_size):
                for y in range(vals_size):
                    vals[y,x] = spline_method(scaled_z,
                                              float(2*y),
                                              float(2*x))
        else:
            for x in range(vals_size):
                for y in range(vals_size):
                    vals[y,x] = spline_method(scaled_z,
                                              float(2*y) + 1.0,
                                              float(2*x) + 1.0)
        return vals

    def getPSFCR(self, z_value):
        return self.getSplineVals(self.spline.f, z_value)

    def getDx(self, z_value):
        return self.getSplineVals(self.spline.dxf, z_value)

    def getDy(self, z_value):
        return self.getSplineVals(self.spline.dyf, z_value)

    def getDz(self, z_value):
        return self.getSplineVals(self.spline.dzf, z_value)    
    

class CRBound3D(object):
    """
    Class for calculating a 3D Cramer-Rao bounds given a spline.
    """
    def __init__(self, spline_file, pixel_size = 160.0):

        self.xy_pixel_size = pixel_size
        self.s_to_psf = CRSplineToPSF3D(pickle.load(open(spline_file, 'rb')))

        self.z_pixel_size = (self.s_to_psf.getZMax() - self.s_to_psf.getZMin())/(float(self.s_to_psf.spline_size) - 1.0)

    def calcCRBound(self, background, photons, z_position = 0.0):
        """
        Calculate Cramer-Rao bounds for a 3D spline.

        Note that this returns the standard deviation 
        bound, not the variance bound.
        """
        # Calculate PSF and it's derivatives.
        psf = self.s_to_psf.getPSFCR(z_position)
        psf_dx = self.s_to_psf.getDx(z_position)
        psf_dy = self.s_to_psf.getDy(z_position)
        psf_dz = self.s_to_psf.getDz(z_position)

        # Normalize to unity & multiply by the number of photons.
        print("psf norm", numpy.sum(psf))
        psf_norm = 1.0/numpy.sum(psf)

        psf_di = psf * psf_norm
        psf_dx = psf_dx * psf_norm * photons / self.xy_pixel_size
        psf_dy = psf_dy * psf_norm * photons / self.xy_pixel_size
        psf_dz = psf_dz * psf_norm * photons / self.z_pixel_size
        psf_dbg = numpy.ones(psf.shape)

        psf_inv = 1.0/(psf_di * photons + background)

        # Calculate Fisher information matrix.
        fmat = numpy.zeros((5,5))
        drvs = [psf_di, psf_dx, psf_dy, psf_dz, psf_dbg]
        for t1 in range(5):
            for t2 in range(5):
                fmat[t1, t2] = numpy.sum(psf_inv * drvs[t1] * drvs[t2])

        fmat_inv = numpy.linalg.inv(fmat)
        crlb = numpy.zeros(5)
        for i in range(5):
            crlb[i] = math.sqrt(fmat_inv[i,i])

            #crlb[1:4] = crlb[1:4] * self.pixel_size
        
        print(crlb)

    def check(self):
        """
        Check that both PSF calculations agree..
        """
        psf_cr = self.s_to_psf.getPSFCR(0.0)
        psf_stp = self.s_to_psf.getPSF(0.0, normalize = False)

        print(numpy.sum(psf_cr), numpy.sum(psf_stp))

    
if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = '(3D) Cramer-Rao bounds calculation')

    parser.add_argument('--spline', dest='spline', type=str, required=True,
                        help = "The name of the spline file")
    parser.add_argument('--background', dest='background', type=float, required=True,
                        help = "The non-specific fluorescence background.")
    parser.add_argument('--photons', dest='photons', type=int, required=True,
                        help = "The number of photons in the PSF image.")

    args = parser.parse_args()
    
    crb = CRBound3D(args.spline)
    crb.check()
    crb.calcCRBound(args.background, args.photons)
