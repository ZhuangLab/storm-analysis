#!/usr/bin/env python
"""
A Python pupil function object. This behaves basically the
same as spliner.spline_to_psf.SplineToPSF object so that
they can be used interchangeably.

Hazen 10/17
"""
import pickle
import numpy

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.simulator.pupil_math as pupilMath

import storm_analysis.pupilfn.pupil_function_c as pupilFnC

class PupilFunction(fitting.PSFFunction):

    def __init__(self, pf_filename = None, zmin = None, zmax = None, **kwds):
        """
        Technically a pupil function would cover any z range, but in 
        fitting we limit it to a finite range. Also, zmin and zmax should 
        be specified in nanometers.
        """
        super(PupilFunction, self).__init__(**kwds)
        self.zmax = zmax
        self.zmin = zmin

        # Load the pupil function data.
        with open(pf_filename, 'rb') as fp:
            pf_data = pickle.load(fp)

        # Get the pupil function and verify that the type is correct.
        pf = pf_data['pf']
        assert (pf.dtype == numpy.complex128)

        # Get the pupil function size and pixel size.
        self.pixel_size = pf_data['pixel_size']
        self.pupil_size = pf.shape[0]

        # Create geometry object.
        if pf_data.get("geo_sim_pf", False):
            geo = pupilMath.GeometrySim(self.pupil_size,
                                        self.pixel_size,
                                        pf_data['wavelength'],
                                        pf_data['immersion_index'],
                                        pf_data['numerical_aperture'])
        else:
            geo = pupilMath.Geometry(self.pupil_size,
                                     self.pixel_size,
                                     pf_data['wavelength'],
                                     pf_data['immersion_index'],
                                     pf_data['numerical_aperture'])

        # Create C pupil function object.
        self.pupil_fn_c = pupilFnC.PupilFunction(geo)
        self.pupil_fn_c.setPF(pf)

    def getCPointer(self):
        return self.pupil_fn_c.getCPointer()

    def getMargin(self):
        return int(self.getSize()/2 + 2)

    def getPSF(self, z_value, shape = None, normalize = False):
        """
        Z value is expected to be in nanometers.
        """
        # Convert z_value to microns.
        z_value = 1.0e-3 * z_value
        
        # Translate to the correct z value.
        self.pupil_fn_c.translate(0.0, 0.0, z_value)

        # Get the (complex) PSF.
        psf = self.pupil_fn_c.getPSF()

        # Convert to intensity
        psf = pupilMath.intensity(psf)

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
        
    def getPixelSize(self):
        return self.pixel_size

    def getScaledZ(self, z):
        return z * 1.0e-3
        
    def getSize(self):
        return self.pupil_size

    
