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
                              ndpointer(dtype = numpy.float64)]

pupil_fn.pfInitialize.argtypes = [ndpointer(dtype = numpy.float64),
                                  ndpointer(dtype = numpy.float64),
                                  ndpointer(dtype = numpy.float64),
                                  ndpointer(dtype = numpy.float64),
                                  ndpointer(dtype = numpy.float64),
                                  ndpointer(dtype = numpy.float64),
                                  ctypes.c_int]
pupil_fn.pfInitialize.restype = ctypes.c_void_p

pupil_fn.pfSetPf.argtypes = [ctypes.c_void_p,
                             ndpointer(dtype = numpy.float64),
                             ndpointer(dtype = numpy.float64)]


class PupilFunction(object):

    def __init__(self, geometry = None, **kwds):
        """
        geometry is a simulation.pupil_math.Geometry object.
        """
        super(PupilFunction, self).__init__(**kwds)
        assert isinstance(geometry, pupilMath.Geometry)

        self.size = geometry.size

        kx = numpy.exp(-1j * 2.0 * numpy.pi * geometry.kx)
        ky = numpy.exp(-1j * 2.0 * numpy.pi * geometry.ky)
        kz = numpy.exp(1j * 2.0 * numpy.pi * geometry.kz)
        
        self.pfn = pupil_fn.pfInitialize(numpy.ascontiguousarray(numpy.real(kx), dtype = numpy.float64),
                                         numpy.ascontiguousarray(numpy.imag(kx), dtype = numpy.float64),
                                         numpy.ascontiguousarray(numpy.real(ky), dtype = numpy.float64),
                                         numpy.ascontiguousarray(numpy.imag(ky), dtype = numpy.float64),
                                         numpy.ascontiguousarray(numpy.real(kz), dtype = numpy.float64),
                                         numpy.ascontiguousarray(numpy.imag(kz), dtype = numpy.float64),
                                         geometry.size)

    def cleanup(self):
        pupil_fn.pfCleanup(self.pfn)
        self.pfn = None
        
    def getPSF(self):
        psf = numpy.zeros((self.size, self.size), dtype = numpy.float64)
        pupil_fn.pfGetPSF(self.pfn, psf)
        return psf
        
    def setPF(self, pf):
        pupil_fn.pfSetPf(self.pfn,
                         numpy.ascontiguousarray(numpy.real(pf), dtype = numpy.float64),
                         numpy.ascontiguousarray(numpy.imag(pf), dtype = numpy.float64))

        
