#!/usr/bin/env python
"""
PSF FFT object for Cramer-Rao bounds calculations.

Hazen 10/17
"""
import numpy
import pickle

import storm_analysis.simulator.pupil_math as pupilMath
import storm_analysis.spliner.cramer_rao as cramerRao

import storm_analysis.pupilfn.pupil_function_c as pupilFnC

class CRPupilFn(cramerRao.CRPSFObject):
    """
    A pupil function based PSF Object for Cramer-Rao bounds calculations.

    The pupil function does not have an intrinsic z range, unlike splines or
    PSF FFTs, so we need to be told what they are.
    """
    def __init__(self, psf_filename = None, zmax = None, zmin = None, **kwds):
        kwds["psf_filename"] = psf_filename
        super(CRPupilFn, self).__init__(**kwds)

        # Load the pupil function data.
        with open(psf_filename, 'rb') as fp:
            pf_data = pickle.load(fp)

        # Get the pupil function and verify that the type is correct.
        pf = pf_data['pf']
        assert (pf.dtype == numpy.complex128)

        # Get the pupil function pixel size.
        self.pupil_size = pf.shape[0]

        # Create geometry object.
        geo = pupilMath.Geometry(self.pupil_size,
                                 pf_data['pixel_size'],
                                 pf_data['wavelength'],
                                 pf_data['immersion_index'],
                                 pf_data['numerical_aperture'])
        
        # Create C pupil function object.
        self.pupil_fn_c = pupilFnC.PupilFunction(geo)
        self.pupil_fn_c.setPF(pf)

        # Additional initializations.
        self.zmax = zmax
        self.zmin = zmin

        # CR weights approximately every 25nm.
        self.n_zvals = int(round((self.zmax - self.zmin)/25.0))
        
        self.delta_xy = self.pixel_size
        self.delta_z = 1000.0

    def cleanup(self):
        self.pupil_fn_c.cleanup()

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
                
    def getDy(self, z_value):
        self.translate(z_value)
        psf_c = self.pupil_fn_c.getPSF()
        psf_c_dx = self.pupil_fn_c.getPSFdx()
        return -2.0 * (numpy.real(psf_c)*numpy.real(psf_c_dx) + numpy.imag(psf_c)*numpy.imag(psf_c_dx))

    def getDx(self, z_value):
        self.translate(z_value)
        psf_c = self.pupil_fn_c.getPSF()
        psf_c_dy = self.pupil_fn_c.getPSFdy()
        return -2.0 * (numpy.real(psf_c)*numpy.real(psf_c_dy) + numpy.imag(psf_c)*numpy.imag(psf_c_dy))
    
    def getDz(self, z_value):
        self.translate(z_value)
        psf_c = self.pupil_fn_c.getPSF()
        psf_c_dz = self.pupil_fn_c.getPSFdz()
        return -2.0 * (numpy.real(psf_c)*numpy.real(psf_c_dz) + numpy.imag(psf_c)*numpy.imag(psf_c_dz))

    def getNormalization(self):
        return self.normalization
        
    def getNZValues(self):
        return self.n_zvals
        
    def getPSF(self, z_value):
        self.translate(z_value)
        return pupilMath.intensity(self.pupil_fn_c.getPSF())
    
    def getZMax(self):
        return self.zmax
    
    def getZMin(self):
        return self.zmin

    def translate(self, z_value):
        self.pupil_fn_c.translate(0.0, 0.0, -z_value * 1.0e-3)

