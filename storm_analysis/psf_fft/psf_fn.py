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


class PSFFn(fitting.PSFFunction):

    def __init__(self, psf_filename = None, **kwds):
        super(PSFFn, self).__init__(**kwds)

        # Load the PSF data.
        with open(psf_filename, 'rb') as fp:
            psf_data = pickle.load(fp)

        # Initialize C library.
        psf = psf_data['psf']
        self.psf_fft_c = psfFFTC.PSFFFT(psf)

        # Store some additional properties.
        self.pixel_size = psf["pixel_size"]
        self.psf_shape = psf.shape
        self.zmax = psf["zmax"]
        self.zmin = psf["zmin"]

        # Sanity checks.
        assert(psf.shape[1] == psf.shape[2])
        assert(self.zmax == -self.zmin)

    def getCPointer(self):
        return self.psf_fft_c.getCPointer()

    def getPixelSize(self):
        return self.pixel_size
        
    def getPSF(self, z_value, shape = None, normalize = False):
        """
        Z value is expected to be in microns.
        """
        # Translate to the correct z value.
        self.psf_fft_c.translate(0.0, 0.0, self.getScaledZ(z_value))

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
        return z_value * (float(self.getZSize()) - 1.0) / (self.zmax - self.zmin)
    
    def getSize(self):
        return self.psf_shape[1]
    
    def getZSize(self):
        return self.psf_shape[0]
        
    def rescaleZ(self):
        return z_value * (self.z_max - self.zmin) / (float(self.getZSize()) - 1.0)
    
