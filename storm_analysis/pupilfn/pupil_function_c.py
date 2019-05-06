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


pupil_fn = loadclib.loadCLibrary("pupil_function")

pupil_fn.pfnCleanup.argtypes = [ctypes.c_void_p]

pupil_fn.pfnGetPSF.argtypes = [ctypes.c_void_p,
                               ndpointer(dtype = numpy.float64),
                               ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPSFIntensity.argtypes = [ctypes.c_void_p,
                                        ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPSFdx.argtypes = [ctypes.c_void_p,
                                 ndpointer(dtype = numpy.float64),
                                 ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPSFdy.argtypes = [ctypes.c_void_p,
                                 ndpointer(dtype = numpy.float64),
                                 ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPSFdz.argtypes = [ctypes.c_void_p,
                                 ndpointer(dtype = numpy.float64),
                                 ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPXEX.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPXEY.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPYEX.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPYEY.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPZEX.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]

pupil_fn.pfnGetPZEY.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]

pupil_fn.pfnInitialize.argtypes = [ndpointer(dtype = numpy.float64),
                                   ndpointer(dtype = numpy.float64),
                                   ndpointer(dtype = numpy.float64),
                                   ctypes.c_int]
pupil_fn.pfnInitialize.restype = ctypes.c_void_p

pupil_fn.pfnSetPF.argtypes = [ctypes.c_void_p,
                              ndpointer(dtype = numpy.float64),
                              ndpointer(dtype = numpy.float64)]

pupil_fn.pfnSetPNEN.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64),
                                ndpointer(dtype = numpy.float64)]
                                
pupil_fn.pfnTranslate.argtypes = [ctypes.c_void_p,
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_double]

pupil_fn.pfnTranslateZ.argtypes = [ctypes.c_void_p,
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

        self.size = geometry.size
        # Size must be an even number.
        assert ((self.size%2)==0)

        # geometry.kz will be a complex number, but the magnitude of the
        # imaginary component is zero so we just ignore it.
        self.pfn = pupil_fn.pfnInitialize(numpy.ascontiguousarray(geometry.kx, dtype = numpy.float64),
                                          numpy.ascontiguousarray(geometry.ky, dtype = numpy.float64),
                                          numpy.ascontiguousarray(numpy.real(geometry.kz), dtype = numpy.float64),
                                          geometry.size)

        if hasattr(geometry, "px_ex"):
            pupil_fn.pfnSetPNEN(self.pfn,
                                numpy.ascontiguousarray(geometry.px_ex, dtype = numpy.float64),
                                numpy.ascontiguousarray(geometry.px_ey, dtype = numpy.float64),
                                numpy.ascontiguousarray(geometry.py_ex, dtype = numpy.float64),
                                numpy.ascontiguousarray(geometry.py_ey, dtype = numpy.float64),
                                numpy.ascontiguousarray(geometry.pz_ex, dtype = numpy.float64),
                                numpy.ascontiguousarray(geometry.pz_ey, dtype = numpy.float64))

    def cleanup(self):
        pupil_fn.pfnCleanup(self.pfn)
        self.pfn = None

    def getCPointer(self):
        return self.pfn
        
    def getPSF(self):
        return self.getXX(pupil_fn.pfnGetPSF)

    def getPSFIntensity(self):
        psf = numpy.zeros((self.size, self.size), dtype = numpy.float64)
        pupil_fn.pfnGetPSFIntensity(self.pfn, psf)
        return psf
    
    def getPSFdx(self):
        return self.getXX(pupil_fn.pfnGetPSFdx)

    def getPSFdy(self):
        return self.getXX(pupil_fn.pfnGetPSFdy)

    def getPSFdz(self):
        return self.getXX(pupil_fn.pfnGetPSFdz)

    def getPXEX(self):
        return self.getXX(pupil_fn.pfnGetPXEX)

    def getPXEY(self):
        return self.getXX(pupil_fn.pfnGetPXEY)
    
    def getPYEX(self):
        return self.getXX(pupil_fn.pfnGetPYEX)
    
    def getPYEY(self):
        return self.getXX(pupil_fn.pfnGetPYEY)

    def getPZEX(self):
        return self.getXX(pupil_fn.pfnGetPZEX)

    def getPZEY(self):
        return self.getXX(pupil_fn.pfnGetPZEY)
    
    def getXX(self, fn):
        r = numpy.zeros((self.size, self.size), dtype = numpy.float64)
        c = numpy.zeros((self.size, self.size), dtype = numpy.float64)
        fn(self.pfn, r, c)
        return r + 1j*c    
    
    def setPF(self, pf):
        pupil_fn.pfnSetPF(self.pfn,
                          numpy.ascontiguousarray(numpy.real(pf), dtype = numpy.float64),
                          numpy.ascontiguousarray(numpy.imag(pf), dtype = numpy.float64))

    def translate(self, dx, dy, dz):
        pupil_fn.pfnTranslate(self.pfn, dx, dy, dz)

    def translateZ(self, dz):
        pupil_fn.pfnTranslateZ(self.pfn, dz)
