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

otf_sc.initialize.argtypes = [ctypes.c_int,
                              ctypes.c_int]

otf_sc.initialize.restype = ctypes.c_void_p

otf_sc.scale.argtypes = [ctypes.c_void_p,
                         ndpointer(dtype = numpy.float64),
                         ndpointer(dtype = numpy.float64)]

otf_sc.setScale.argtypes = [ctypes.c_void_p,
                            ndpointer(dtype = numpy.float64)]


class OTFScaler(object):

    def __init__(self, fftw_estimate = False, size = None, **kwds):
        """
        fftw_estimate - FFTW estimates best way to perform FFT.
        size - Size in pixels of the PSF to scale.
        """
        self.size = size
        self.scaler = otf_sc.initialize(self.size, int(fftw_estimate))

    def cleanup(self):
        otf_sc.cleanup(self.scaler)

    def scale(self, psf):
        result = numpy.zeros(psf.shape, dtype = numpy.float64)
        otf_sc.scale(self.scaler,
                     numpy.ascontiguousarray(psf, dtype = numpy.float64),
                     result)
        return result

    def setScale(self, scale):
        """
        scale - An array of numpy floats, symmetric in X/Y.
        """
        assert (scale.shape[0] == self.size)
        assert (scale.shape[1] == self.size)

        scale = numpy.fft.fftshift(scale)
        otf_sc.setScale(self.scaler, numpy.ascontiguousarray(scale, dtype = numpy.float64))
        
