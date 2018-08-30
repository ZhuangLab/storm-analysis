#!/usr/bin/env python
"""
Simple Python interface to psf_fft.c

Hazen 10/17
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

import storm_analysis.sa_library.loadclib as loadclib

psf_fft = loadclib.loadCLibrary("psf_fft")

psf_fft.pFTCleanup.argtypes = [ctypes.c_void_p]

psf_fft.pFTGetPSF.argtypes = [ctypes.c_void_p,
                              ndpointer(dtype = numpy.float64)]

psf_fft.pFTGetPSFdx.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64)]

psf_fft.pFTGetPSFdy.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64)]

psf_fft.pFTGetPSFdz.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64)]

psf_fft.pFTInitialize.argtypes = [ndpointer(dtype = numpy.float64),
                                  ctypes.c_int,
                                  ctypes.c_int,
                                  ctypes.c_int]
psf_fft.pFTInitialize.restype = ctypes.c_void_p

psf_fft.pFTTranslate.argtypes = [ctypes.c_void_p,
                                 ctypes.c_double,
                                 ctypes.c_double,
                                 ctypes.c_double]


class PSFFFTException(Exception):

    def __init__(self, message):
        Exception.__init__(self, message)

        
class PSFFFT(object):

    def __init__(self, psf = None, **kwds):
        super(PSFFFT, self).__init__(**kwds)
        
        self.psf_shape = psf.shape

        c_psf = numpy.ascontiguousarray(psf, dtype = numpy.float64)
        self.pfft = psf_fft.pFTInitialize(c_psf,
                                          self.psf_shape[0],
                                          self.psf_shape[1],
                                          self.psf_shape[2])

    def cleanup(self):
        psf_fft.pFTCleanup(self.pfft)

    def getCPointer(self):
        return self.pfft
    
    def getPSF(self):
        psf = numpy.ascontiguousarray(numpy.zeros((self.psf_shape[1], self.psf_shape[2]), dtype = numpy.float64),
                                      dtype = numpy.float64)
        psf_fft.pFTGetPSF(self.pfft, psf)
        return psf

    def getPSFdx(self):
        dx = numpy.ascontiguousarray(numpy.zeros((self.psf_shape[1], self.psf_shape[2]), dtype = numpy.float64),
                                      dtype = numpy.float64)
        psf_fft.pFTGetPSFdx(self.pfft, dx)
        return dx

    def getPSFdy(self):
        dy = numpy.ascontiguousarray(numpy.zeros((self.psf_shape[1], self.psf_shape[2]), dtype = numpy.float64),
                                      dtype = numpy.float64)
        psf_fft.pFTGetPSFdy(self.pfft, dy)
        return dy

    def getPSFdz(self):
        dz = numpy.ascontiguousarray(numpy.zeros((self.psf_shape[1], self.psf_shape[2]), dtype = numpy.float64),
                                      dtype = numpy.float64)
        psf_fft.pFTGetPSFdz(self.pfft, dz)
        return dz

    def translate(self, dx, dy, dz):
        """
        Note:
          1. dx and dy are reversed and have the opposite sign from simulator.pupil_math.
          2. dz is in z pixel units, not microns.
        """
        psf_fft.pFTTranslate(self.pfft, dx, dy, dz)
        
