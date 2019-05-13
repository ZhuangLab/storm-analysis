#!/usr/bin/env python
"""
A Python PSF FFT function object.

Note: 
  1. This will not work directly with the measured PSFs created by
     Spliner and/or Multifit, as these are 2x upsampled. You need
     to downsample them first with ..

  2. The range covered by the measured PSF must be symmetric in Z,
     i.e. it goes from -z_range to z_range.

Hazen 10/17
"""
import pickle
import numpy

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.simulator.pupil_math as pupilMath

import storm_analysis.psf_fft.psf_fft_c as psfFFTC


class PSFFnBase(fitting.PSFFunction):

    def __init__(self, psf_data = None, **kwds):
        super(PSFFnBase, self).__init__(**kwds)

        # Initialize C library.
        psf = psf_data['psf']
        self.psf_fft_c = psfFFTC.PSFFFT(psf)

        # Store some additional properties.
        self.pixel_size = psf_data["pixel_size"]
        self.psf_shape = psf.shape

        # These are in units of nanometers.
        self.zmax = psf_data["zmax"] 
        self.zmin = psf_data["zmin"]

        self.scale_gSZ = (float(self.getZSize()) - 1.0) / (self.zmax - self.zmin)
        self.scale_rZ = 1.0e-3 * (self.zmax - self.zmin) / (float(self.getZSize()) - 1.0)
        
        # Sanity checks.
        assert ((psf.shape[0]%2) == 1), "Z size must be an odd number."
        assert (psf.shape[1] == psf.shape[2]), "X/Y size must be the same."
        assert (self.zmax == -self.zmin), "z range must be symmetric."

    def getCPointer(self):
        return self.psf_fft_c.getCPointer()

    def getMargin(self):
        return int((self.getSize() + 1)/2 + 2)

    def getPixelSize(self):
        return self.pixel_size
        
    def getPSF(self, z_value, shape = None, normalize = False):
        """
        Z value is expected to be in nanometers.
        """
        # Translate to the correct x/y/z value.
        #
        # Why 1.0, 1.0 offset in X/Y? We do this so that the PSF will match
        # that of Spliner (spline_to_psf.SplineToPSF3D.getPSF()).
        #
        self.psf_fft_c.translate(1.0, 1.0, self.getScaledZ(z_value))

        # Get the (complex) PSF.
        psf = self.psf_fft_c.getPSF()

        # Center into a (larger) array if requested.
        if shape is not None:
            psf_size = psf.shape[0]
            im_size_x = shape[0]
            im_size_y = shape[1]

            start_x = int(im_size_x/2.0 - psf_size/2.0)
            start_y = int(im_size_y/2.0 - psf_size/2.0)

            end_x = start_x + psf_size
            end_y = start_y + psf_size

            temp = numpy.zeros((im_size_x, im_size_y))
            temp[start_x:end_x,start_y:end_y] = psf

            psf = temp

        # Normalize if requested.
        if normalize:
            psf = psf/numpy.sum(psf)

        return psf

    def getScaledZ(self, z_value):
        """
        This expects z_value to be in nanometers.
        """
        return z_value * self.scale_gSZ
    
    def getSize(self):
        return self.psf_shape[1]
    
    def getZSize(self):
        return self.psf_shape[0]
        
    def rescaleZ(self, z_value):
        """
        This returns a z_value in microns.
        """
        return z_value * self.scale_rZ


class PSFFn(PSFFnBase):

    def __init__(self, psf_filename = None, **kwds):

        # Load the PSF data.
        with open(psf_filename, 'rb') as fp:
            psf_data = pickle.load(fp)

        super(PSFFn, self).__init__(psf_data = psf_data, **kwds)
        
