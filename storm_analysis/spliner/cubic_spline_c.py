#!/usr/bin/env python
"""
Simple Python interface to cubic_spline.c. This is used mostly for 
testing that the C library works correctly.

Hazen 12/13
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import random
import sys

import storm_analysis.spliner.spline2D as spline2D
import storm_analysis.spliner.spline3D as spline3D

import storm_analysis.sa_library.loadclib as loadclib

# Load the library.
cubic = loadclib.loadCLibrary("cubic_spline")

# C interface definition.
cubic.computeDelta2D.argtypes = [ctypes.c_void_p,
                                 ctypes.c_double,
                                 ctypes.c_double]

cubic.computeDelta3D.argtypes = [ctypes.c_void_p,
                                 ctypes.c_double,
                                 ctypes.c_double,
                                 ctypes.c_double]

cubic.dxfSpline2D.argtypes = [ctypes.c_void_p,
                              ctypes.c_double,
                              ctypes.c_double]
cubic.dxfSpline2D.restype = ctypes.c_double

cubic.dxfSpline3D.argtypes = [ctypes.c_void_p,
                              ctypes.c_double,
                              ctypes.c_double,
                              ctypes.c_double]
cubic.dxfSpline3D.restype = ctypes.c_double

cubic.dyfSpline2D.argtypes = [ctypes.c_void_p,
                              ctypes.c_double,
                              ctypes.c_double]
cubic.dyfSpline2D.restype = ctypes.c_double

cubic.dyfSpline3D.argtypes = [ctypes.c_void_p,
                              ctypes.c_double,
                              ctypes.c_double,
                              ctypes.c_double]
cubic.dyfSpline3D.restype = ctypes.c_double

cubic.dzfSpline3D.argtypes = [ctypes.c_void_p,
                              ctypes.c_double,
                              ctypes.c_double,
                              ctypes.c_double]
cubic.dzfSpline3D.restype = ctypes.c_double

cubic.fSpline2D.argtypes = [ctypes.c_void_p,
                            ctypes.c_double,
                            ctypes.c_double]
cubic.fSpline2D.restype = ctypes.c_double

cubic.fSpline3D.argtypes = [ctypes.c_void_p,
                            ctypes.c_double,
                            ctypes.c_double,
                            ctypes.c_double]
cubic.fSpline3D.restype = ctypes.c_double

cubic.getPSF2D.argtypes = [ctypes.c_void_p,
                           ndpointer(dtype=numpy.float64),
                           ctypes.c_double,
                           ctypes.c_double]

cubic.getPSF3D.argtypes = [ctypes.c_void_p,
                           ndpointer(dtype=numpy.float64),
                           ctypes.c_double,
                           ctypes.c_double,
                           ctypes.c_double]

cubic.initSpline2D.argtypes = [ndpointer(dtype=numpy.float64),
                               ctypes.c_int,
                               ctypes.c_int]
cubic.initSpline2D.restype = ctypes.c_void_p

cubic.initSpline3D.argtypes = [ndpointer(dtype=numpy.float64),
                               ctypes.c_int,
                               ctypes.c_int,
                               ctypes.c_int]
cubic.initSpline3D.restype = ctypes.c_void_p

cubic.splineCleanup.argtypes = [ctypes.c_void_p]


class CubicSplineCException(Exception):
    pass


# Classes.

class CSpline(object):

    def checkCSpline(self):
        if self.c_spline is None:
            raise CubicSplineCException("Pointer to C spline data structure is NULL")

    def cleanup(self):
        self.checkCSpline()
        cubic.splineCleanup(self.c_spline)
        self.c_spline = None

    def getCPointer(self):
        return self.c_spline
        
        
class CSpline2D(CSpline):

    def __init__(self, py_spline):
        self.py_spline = py_spline
        self.c_spline = cubic.initSpline2D(numpy.ascontiguousarray(self.py_spline.coeff, dtype = numpy.float64),
                                           self.py_spline.max_i,
                                           self.py_spline.max_i)

    def dxf(self, y, x):
        self.checkCSpline()
        return cubic.dxfSpline2D(self.c_spline, y, x)

    def dyf(self, y, x):
        self.checkCSpline()
        return cubic.dyfSpline2D(self.c_spline, y, x)

    def f(self, y, x):
        self.checkCSpline()
        return cubic.fSpline2D(self.c_spline, y, x)

    def getPSF(self, dy, dx):
        psf_size = self.py_spline.max_i - 1
        psf = numpy.zeros((psf_size, psf_size), dtype = numpy.float64)
        cubic.getPSF2D(self.c_spline, psf, dy, dx)
        return psf

    def py_f(self, y, x):
        return self.py_spline.f(x, y)


class CSpline3D(CSpline):

    def __init__(self, py_spline):
        self.py_spline = py_spline
        self.c_spline = cubic.initSpline3D(numpy.ascontiguousarray(self.py_spline.coeff, dtype = numpy.float64),
                                           self.py_spline.max_i,
                                           self.py_spline.max_i,
                                           self.py_spline.max_i)
        
    def dxf(self, z, y, x):
        self.checkCSpline() 
        return cubic.dxfSpline3D(self.c_spline, z, y, x)

    def dyf(self, z, y, x):
        self.checkCSpline()
        return cubic.dyfSpline3D(self.c_spline, z, y, x)
        
    def dzf(self, z, y, x):
        self.checkCSpline()
        return cubic.dzfSpline3D(self.c_spline, z, y, x)
        
    def f(self, z, y, x):
        self.checkCSpline()
        return cubic.fSpline3D(self.c_spline, z, y, x)

    def getPSF(self, z, dy, dx):
        psf_size = self.py_spline.max_i - 1
        psf = numpy.zeros((psf_size, psf_size), dtype = numpy.float64)
        cubic.getPSF3D(self.c_spline, psf, z, dy, dx)
        return psf
        
    def py_f(self, z, y, x):
        return self.py_spline.f(z, y, x)


