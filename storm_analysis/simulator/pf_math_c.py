#!/usr/bin/env python
"""
Simple Python interface to zernike.c

Hazen 10/14
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.loadclib as loadclib

pf_math = loadclib.loadCLibrary("pf_math")

pf_math.pfCleanup.argtypes = [ctypes.c_void_p]

pf_math.pfGetPSF.argtypes = [ctypes.c_void_p,
                             ndpointer(dtype=numpy.float64),
                             ctypes.c_double]

pf_math.pfInitialize.argtypes = [ndpointer(dtype=numpy.float64),
                                 ctypes.c_int]
pf_math.pfInitialize.restype = ctypes.c_void_p

pf_math.pfResetScaling.argtypes = [ctypes.c_void_p]

pf_math.pfSetPF.argtypes = [ctypes.c_void_p,
                            ndpointer(dtype=numpy.float64),
                            ndpointer(dtype=numpy.float64)]

pf_math.pfSetScaling.argtypes = [ctypes.c_void_p,
                                 ndpointer(dtype=numpy.float64)]
                            
pf_math.pfZernike.argtypes = [ctypes.c_int,
                              ctypes.c_int,
                              ctypes.c_double,
                              ctypes.c_double]
pf_math.pfZernike.restype = ctypes.c_double

pf_math.pfZernikeGrid.argtypes = [ndpointer(dtype=numpy.float64),
                                  ctypes.c_int,
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_int,
                                  ctypes.c_int]


class PupilMath(object):
    """
    This speeds up the calculation of the PSF from the PF, which is what
    takes 90-95% of the time in pupil function calculations.
    """
    def __init__(self, kz = None, **kwds):
        super(PupilMath, self).__init__(**kwds)

        assert (len(kz.shape) == 2), "kz must be two dimensional."
        assert (kz.shape[0] == kz.shape[1]), "kz must be square."

        self.pf_size = kz.shape[0]
        assert ((self.pf_size%2) == 0), "pf size must be an even number."
        
        c_kz = numpy.ascontiguousarray(numpy.real(kz), dtype = numpy.float64)
        self.pfm = pf_math.pfInitialize(c_kz, kz.shape[0])

    def cleanup(self):
        pf_math.pfCleanup(self.pfm)
        self.pfm = None

    def getPSF(self, z_vals):
        psf = numpy.ascontiguousarray(numpy.zeros((len(z_vals), self.pf_size, self.pf_size)),
                                      dtype = numpy.float64)

        for i, z in enumerate(z_vals):
            pf_math.pfGetPSF(self.pfm, psf[i,:,:], z)

        return psf
               
    def setPF(self, pf):
        assert(pf.shape[0] == self.pf_size), "PF has the wrong shape, axis 0."
        assert(pf.shape[1] == self.pf_size), "PF has the wrong shape, axis 1."

        #pf = numpy.fft.fftshift(pf)
        pf_r = numpy.ascontiguousarray(numpy.real(pf), dtype = numpy.float64)
        pf_c = numpy.ascontiguousarray(numpy.imag(pf), dtype = numpy.float64)
        pf_math.pfSetPF(self.pfm, pf_r, pf_c)

    def setScaling(self, scaling):
        if scaling is None:
            pf_math.pfResetScaling(self.pfm)

        else:
            assert(scaling.shape[0] == self.pf_size), "Scaling has the wrong shape, axis 0."
            assert(scaling.shape[1] == self.pf_size), "Scaling has the wrong shape, axis 1."

            scaling = numpy.fft.fftshift(scaling)
            c_scaling = numpy.ascontiguousarray(scaling, dtype = numpy.float64)
            pf_math.pfSetScaling(self.pfm, scaling)


def zernikeGrid(np_array, scale, m, n, radius = None, center = None):
    """
    Calculate the shape of a Zernike polynomial and store it
    in np_array.
    """
    if (np_array.shape[0] != np_array.shape[1]):
        print("Array must be square.")
        return
    
    if radius is None:
        radius = np_array.shape[0]/2
        
    if center is None:
        #center = (np_array.shape[0]/2 - 0.5)
        center = np_array.shape[0] * 0.5

    c_np_array = numpy.ascontiguousarray(np_array, dtype=numpy.float64)

    pf_math.pfZernikeGrid(c_np_array,
                          np_array.shape[0],
                          center,
                          radius,
                          scale,
                          m,
                          n)
    
    return c_np_array
