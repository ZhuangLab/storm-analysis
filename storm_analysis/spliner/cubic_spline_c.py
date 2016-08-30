#!/usr/bin/python
#
# Simple Python interface to cubic_spline.c. This is used mostly for 
# testing that the C library works correctly.
#
# Hazen 12/13
#

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import random
import sys

import spline2D
import spline3D

import storm_analysis.sa_library.loadclib as loadclib

# Load the library.
cubic = loadclib.loadCLibrary(os.path.dirname(__file__), "cubic_spline")

# C interface definition.
cubic.computeDelta2D.argtypes = [ctypes.c_double,
                                 ctypes.c_double]
cubic.computeDelta3D.argtypes = [ctypes.c_double,
                                 ctypes.c_double,
                                 ctypes.c_double]

#cubic.fAt2D.argtypes = [ctypes.c_int,
#                        ctypes.c_int]
#cubic.fAt2D.restype = ctypes.c_double
#cubic.fAt3D.argtypes = [ctypes.c_int,
#                        ctypes.c_int,
#                        ctypes.c_int]
#cubic.fAt3D.restype = ctypes.c_double
#
cubic.dxfSpline2D.argtypes = [ctypes.c_double,
                              ctypes.c_double]
cubic.dxfSpline2D.restype = ctypes.c_double
cubic.dxfSpline3D.argtypes = [ctypes.c_double,
                              ctypes.c_double,
                              ctypes.c_double]
cubic.dxfSpline3D.restype = ctypes.c_double

cubic.dyfSpline2D.argtypes = [ctypes.c_double,
                              ctypes.c_double]
cubic.dyfSpline2D.restype = ctypes.c_double
cubic.dyfSpline3D.argtypes = [ctypes.c_double,
                              ctypes.c_double,
                              ctypes.c_double]
cubic.dyfSpline3D.restype = ctypes.c_double

cubic.dzfSpline3D.argtypes = [ctypes.c_double,
                              ctypes.c_double,
                              ctypes.c_double]
cubic.dzfSpline3D.restype = ctypes.c_double

cubic.fSpline2D.argtypes = [ctypes.c_double,
                            ctypes.c_double]
cubic.fSpline2D.restype = ctypes.c_double

cubic.fSpline3D.argtypes = [ctypes.c_double,
                            ctypes.c_double,
                            ctypes.c_double]
cubic.fSpline3D.restype = ctypes.c_double

cubic.initSpline2D.argtypes = [ndpointer(dtype=numpy.float64),
                               ctypes.c_int,
                               ctypes.c_int]
cubic.initSpline3D.argtypes = [ndpointer(dtype=numpy.float64),
                               ctypes.c_int,
                               ctypes.c_int,
                               ctypes.c_int]


# Classes.
class CSpline2D():

    def __init__(self, d):

        self.py_spline = spline2D.Spline2D(d)
        cubic.initSpline2D(numpy.ascontiguousarray(self.py_spline.coeff, dtype = numpy.float64),
                           self.py_spline.max_i,
                           self.py_spline.max_i)

    def dxf(self, x, y):
        return cubic.dxfSpline2D(x, y)

    def dyf(self, x, y):
        return cubic.dyfSpline2D(x, y)

    def f(self, x, y):
        return cubic.fSpline2D(x, y)

    def py_f(self, x, y):
        return self.py_spline.f(x, y)

class CSpline3D():

    def __init__(self, d):

        self.py_spline = spline3D.Spline3D(d)
        cubic.initSpline3D(numpy.ascontiguousarray(self.py_spline.coeff, dtype = numpy.float64),
                           self.py_spline.max_i,
                           self.py_spline.max_i,
                           self.py_spline.max_i)

    def dxf(self, x, y, z):
        return cubic.dxfSpline3D(x, y, z)

    def dyf(self, x, y, z):
        return cubic.dyfSpline3D(x, y, z)

    def dzf(self, x, y, z):
        return cubic.dzfSpline3D(x, y, z)

    def f(self, x, y, z):
        return cubic.fSpline3D(x, y, z)

    def py_f(self, x, y, z):
        return self.py_spline.f(x, y, z)


# Tests.
if __name__ == "__main__":

    if 0:
        x = numpy.arange(0.0, 2.001, 0.5 - 1.0e-12)
        y = numpy.array([[0.0, 1.0, 2.0],
                         [3.0, 4.0, 5.0],
                         [6.0, 7.0, 8.0]])

        s = CSpline2D(y)

        if 0:
            for i in range(10):
                px = 2.0 * random.random()
                py = 2.0 * random.random()
                print(s.f(px, py), s.py_f(px, py))
                
        if 1:
            surf = numpy.zeros((x.size, x.size))
            dx_surf = numpy.zeros((x.size, x.size))
            dy_surf = numpy.zeros((x.size, x.size))
            for i in range(x.size):
                for j in range(x.size):
                    surf[i,j] = s.f(x[i],x[j])
                    dx_surf[i,j] = s.dxf(x[i],x[j])
                    dy_surf[i,j] = s.dyf(x[i],x[j])
                
            print(surf)
            print(dx_surf)
            print(dy_surf)

    if 1:
        x = numpy.arange(0.0, 2.01, 1.0 - 1.0e-12)
        y = numpy.array([[[0.0, 1.0, 2.0],
                          [3.0, 4.0, 5.0],
                          [6.0, 7.0, 8.0]],
                         [[9.0, 10.0, 11.0],
                          [12.0, 13.0, 14.0],
                          [15.0, 16.0, 17.0]],
                         [[18.0, 19.0, 20.0],
                          [21.0, 22.0, 23.0],
                          [24.0, 25.0, 26.0]]])

        s = CSpline3D(y)

        if 1:
            surf = numpy.zeros((x.size, x.size, x.size))
            dx_surf = numpy.zeros((x.size, x.size, x.size))
            dy_surf = numpy.zeros((x.size, x.size, x.size))
            dz_surf = numpy.zeros((x.size, x.size, x.size))
            for i in range(x.size):
                for j in range(x.size):
                    for k in range(x.size):
                        surf[i,j,k] = s.f(x[i],x[j],x[k])
                        dx_surf[i,j,k] = s.dxf(x[i],x[j],x[k])
                        dy_surf[i,j,k] = s.dyf(x[i],x[j],x[k])
                        dz_surf[i,j,k] = s.dzf(x[i],x[j],x[k])

            print("f:")
            print(surf)
#            print "dxf:"
#            print dx_surf
#            print "dyf:"
#            print dy_surf
#            print "dzf:"
#            print dz_surf
