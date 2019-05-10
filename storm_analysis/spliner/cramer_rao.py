#!/usr/bin/env python
"""
Estimate Cramer-Rao bounds for a spline based on the formulation in
this paper (specifically equation 18):

"Three dimensional single molecule localization using a phase retrieved pupil
function", Sheng Liu, Emil B. Kromann, Wesley D. Krueger, Joerg Bewersdorf, 
and Keith A. Lidke, Optics Express 2013.

Note: Tested against Pupil Function and PSF FFT with 
      storm_analysis/diagnostics/cramer_rao/

Hazen 02/17
"""
import numpy
import pickle

import storm_analysis.spliner.spline3D as spline3D

# FIXME: Also handle 2D splines.

class CramerRaoException(Exception):
    pass

class CRPSFObject(object):
    """
    Base class for PSF objects for Cramer-Rao bounds calculations.
    """
    def __init__(self, psf_data = None, pixel_size = None, **kwds):
        super(CRPSFObject, self).__init__(**kwds)
        self.pixel_size = pixel_size
        
        if "normalization" in psf_data:
            self.normalization = psf_data["normalization"]
        else:
            self.normalization = 1.0
            print("No normalization data found, using 1.0.")

    def cleanup(self):
        """
        This is useful for C library based implementations.
        """
        pass
    

class CRSplineToPSF3D(CRPSFObject):
    """
    A Spline based PSF Object for Cramer-Rao bounds calculations.
    """
    def __init__(self, psf_filename = None, **kwds):

        # Load the spline.
        with open(psf_filename, 'rb') as fp:
            spline_data = pickle.load(fp)

        super(CRSplineToPSF3D, self).__init__(psf_data = spline_data, **kwds)
            
        self.zmax = spline_data["zmax"]
        self.zmin = spline_data["zmin"]
        self.spline = spline3D.Spline3D(spline_data["spline"], spline_data["coeff"])

        self.delta_xy = self.pixel_size
        self.delta_z = (self.getZMax() - self.getZMin())/float(self.spline.getSize())

    def getDeltaXY(self):
        """
        Return delta XY scaling term (in nanometers).
        """
        return self.delta_xy
        
    def getDeltaZ(self):
        """
        Return delta Z scaling term (in nanometers).
        """
        return self.delta_z
                
    def getDx(self, z_value):
        return self.getSplineVals(self.spline.dxf, z_value)

    def getDy(self, z_value):
        return self.getSplineVals(self.spline.dyf, z_value)

    def getDz(self, z_value):
        return self.getSplineVals(self.spline.dzf, z_value)

    def getNormalization(self):
        return self.normalization
        
    def getNZValues(self):
        return self.spline.getSize() + 1
        
    def getPSF(self, z_value):
        return self.getSplineVals(self.spline.f, z_value)

    def getSplineVals(self, spline_method, z_value):
        """
        Return the spline values of a given method such as:

        spline.f - The PSF
        spline.dx - The derivative in x
        ...

        This is basically a specialized version of the
        spline_to_psf.SplineToPSF3D.getPSF() method.
        """
        scaled_z = float(self.spline.getSize()) * (z_value - self.zmin) / (self.zmax - self.zmin)
        
        vals_size = self.spline.getSize() - 1

        vals = numpy.zeros((vals_size, vals_size))
        if((vals_size%2) == 0):
            for x in range(vals_size):
                for y in range(vals_size):
                    vals[y,x] = spline_method(scaled_z,
                                              float(y) + 0.5,
                                              float(x) + 0.5)
        else:
            for x in range(vals_size):
                for y in range(vals_size):
                    vals[y,x] = spline_method(scaled_z,
                                              float(2*y) + 1.0,
                                              float(2*x) + 1.0)
        return vals

    def getZMax(self):
        #
        # We subtract a small amount so that we don't get complaints
        # about the z value being out of range for the spline when
        # it is exactly at the maximum.
        #
        return self.zmax - 1.0e-9
    
    def getZMin(self):
        return self.zmin + 1.0e-9

    
class CRBound3D(object):
    """
    Class for calculating a 3D Cramer-Rao bounds given a spline.

    Notes: 
      (1) This returns the variance.
      (2) The results for x,y and z are nanometers.

    FIXME: Check that we have not broken the scripts provided in the
           C-Spline paper, which I think is the only reason why we
           are maintaining this class.
    """
    def __init__(self, spline_file, pixel_size = 160.0):
        self.cr_psf_object = CRSplineToPSF3D(spline_file = spline_file,
                                             pixel_size = pixel_size)

    def calcCRBound(self, background, photons, z_position = 0.0):
        """
        This just calls the calcCRBound3D() function.
        """
        return calcCRBound3D(self.cr_psf_object,
                             background,
                             photons,
                             z_position)
    
            
def calcCRBound3D(cr_psf_object, background, photons, z_value):
    """
    Class for calculating a 3D Cramer-Rao bounds given a CR PSF object.

    Notes: 
     (1) This returns the variance.
     (2) The results for x,y and z are nanometers.    
     (3) z_value is expected to be in microns.
    """
    # Convert z_value to nanometers.
    z_value = z_value * 1.0e+3
    
    # Calculate PSF and it's derivatives.
    psf = cr_psf_object.getPSF(z_value)
    psf_dx = cr_psf_object.getDx(z_value)
    psf_dy = cr_psf_object.getDy(z_value)
    psf_dz = cr_psf_object.getDz(z_value)

    # Normalize to unity & multiply by normalization constant. For a single
    # PSF getNormalization() returns 1.0, but for PSFs that have been processed
    # with multiplane/normalize_psfs it may be less than 1.0.
    #
    psf_norm = cr_psf_object.getNormalization()/numpy.sum(psf)

    psf_di = psf * psf_norm        
    psf_dx = -psf_dx * psf_norm * photons / cr_psf_object.getDeltaXY()
    psf_dy = -psf_dy * psf_norm * photons / cr_psf_object.getDeltaXY()
    psf_dz = psf_dz * psf_norm * photons / cr_psf_object.getDeltaZ()
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
        crlb[i] = fmat_inv[i,i]
        
    return crlb

    
if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = '(3D) Cramer-Rao bounds calculation, results in nanometers (for X/Y/Z)')

    parser.add_argument('--spline', dest='spline', type=str, required=True,
                        help = "The name of the spline file")
    parser.add_argument('--background', dest='background', type=float, required=True,
                        help = "The non-specific fluorescence background.")
    parser.add_argument('--photons', dest='photons', type=int, required=True,
                        help = "The number of photons in the PSF image.")
    parser.add_argument('--pixel_size', dest='pixel_size', type=float, required=False, default=160.0,
                        help = "The XY pixel size in nanometers.")

    args = parser.parse_args()

    cr_po = CRSplineToPSF3D(psf_filename = args.spline,
                            pixel_size = args.pixel_size)
    print(calcCRBound3D(cr_po, args.background, args.photons, 0.0))
