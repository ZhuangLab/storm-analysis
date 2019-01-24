#!/usr/bin/env python
"""
Simple Python interface to otf_scaling.c

Hazen 01/19
"""
import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.recenter_psf as recenterPSF


otf_sc = loadclib.loadCLibrary("otf_scaling")

otf_sc.cleanup.argtypes = [ctypes.c_void_p]

otf_sc.initialize.argtypes = [ndpointer(dtype = numpy.float64),
                              ctypes.c_int,
                              ctypes.c_int]

otf_sc.initialize.restype = ctypes.c_void_p

otf_sc.scale.argtypes = [ctypes.c_void_p,
                         ndpointer(dtype = numpy.float64),
                         ndpointer(dtype = numpy.float64)]


class OTFScaler(object):

    def __init__(self, estimate_fft_plan = False, geometry = None, sigma = None, **kwds):
        """
        geometry - A simulation.pupil_math.Geometry object.
        sigma - The sigma of the Gaussian to use for OTF scaling.
        """
        otf_scaler = numpy.fft.fftshift(geometry.gaussianScalingFactor(sigma))
        self.scaler = otf_sc.initialize(numpy.ascontiguousarray(otf_scaler, dtype = numpy.float64),
                                        otf_scaler.shape[0],
                                        int(estimate_fft_plan))

    def cleanup(self):
        otf_sc.cleanup(self.scaler)

    def scale(self, psf):
        result = numpy.zeros(psf.shape, dtype = numpy.float64)
        otf_sc.scale(self.scaler,
                     numpy.ascontiguousarray(psf, dtype = numpy.float64),
                     result)
        return result
    
