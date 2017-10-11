#!/usr/bin/env python
"""
Simple Python interface to pupil_function.c

Hazen 10/17
"""
import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

import storm_analysis.sa_library.loadclib as loadclib

import storm_analysis.simulator.pupil_math as pupilMath

#import storm_analysis.sa_library.recenter_psf as recenterPSF

pupil_fn = loadclib.loadCLibrary("storm_analysis.sa_library", "pupil_function")

pupil_fn.pfCleanup.argtypes = [ctypes.c_void_p]

pupil_fn.pfGetPSF.argtypes = [ctypes.c_void_p,
                              ndpointer(dtype = numpy.float64),
                              ndpointer(dtype = numpy.float64)]

pupil_fn.pfGetPSFdx.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]

pupil_fn.pfGetPSFdy.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]

pupil_fn.pfGetPSFdz.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]

pupil_fn.pfInitialize.argtypes = [ndpointer(dtype = numpy.float64),
                                  ndpointer(dtype = numpy.float64),
                                  ndpointer(dtype = numpy.float64),
                                  ctypes.c_int]
pupil_fn.pfInitialize.restype = ctypes.c_void_p

pupil_fn.pfSetPF.argtypes = [ctypes.c_void_p,
                             ndpointer(dtype = numpy.float64),
                             ndpointer(dtype = numpy.float64)]

pupil_fn.pfTranslate.argtypes = [ctypes.c_void_p,
                                 ctypes.c_double,
                                 ctypes.c_double,
                                 ctypes.c_double]


class PupilFunction(object):

    def __init__(self, geometry = None, **kwds):
        """
        geometry is a simulation.pupil_math.Geometry object.

        Note: All of the getXY functions return complex numbers. In order
              to convert the output of getPSF() to a magnitude you need 
              to do mag = psf*numpy.conj(psf).
        """
        super(PupilFunction, self).__init__(**kwds)
        assert isinstance(geometry, pupilMath.Geometry)

        self.size = geometry.size
        # Size must be an even number.
        assert ((self.size%2)==0)

        # geometry.kz will be a complex number, but the magnitude of the
        # imaginary component is zero so we just ignore it.
        self.pfn = pupil_fn.pfInitialize(numpy.ascontiguousarray(geometry.kx, dtype = numpy.float64),
                                         numpy.ascontiguousarray(geometry.ky, dtype = numpy.float64),
                                         numpy.ascontiguousarray(numpy.real(geometry.kz), dtype = numpy.float64),
                                         geometry.size)

    def cleanup(self):
        pupil_fn.pfCleanup(self.pfn)
        self.pfn = None

    def getPSF(self):
        return self.getXX(pupil_fn.pfGetPSF)
    
    def getPSFdx(self):
        return self.getXX(pupil_fn.pfGetPSFdx)

    def getPSFdy(self):
        return self.getXX(pupil_fn.pfGetPSFdy)

    def getPSFdz(self):
        return self.getXX(pupil_fn.pfGetPSFdz)

    def getXX(self, fn):
        r = numpy.zeros((self.size, self.size), dtype = numpy.float64)
        c = numpy.zeros((self.size, self.size), dtype = numpy.float64)
        fn(self.pfn, r, c)
        return r + 1j*c    
    
    def setPF(self, pf):
        pupil_fn.pfSetPF(self.pfn,
                         numpy.ascontiguousarray(numpy.real(pf), dtype = numpy.float64),
                         numpy.ascontiguousarray(numpy.imag(pf), dtype = numpy.float64))

    def translate(self, dx, dy, dz):
        pupil_fn.pfTranslate(self.pfn, dx, dy, dz)
