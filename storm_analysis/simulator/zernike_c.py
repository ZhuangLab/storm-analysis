#!/usr/bin/python
#
# Simple Python interface to zernike.c
#
# Hazen 10/14
#

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.loadclib as loadclib

zernike = loadclib.loadCLibrary("storm_analysis.simulator", "zernike")

zernike.zernike.argtypes = [ctypes.c_int,
                            ctypes.c_int,
                            ctypes.c_double,
                            ctypes.c_double]
zernike.zernike.restype = ctypes.c_double

zernike.zernike_grid.argtypes = [ndpointer(dtype=numpy.float64),
                                 ctypes.c_int,
                                 ctypes.c_double,
                                 ctypes.c_double,
                                 ctypes.c_double,
                                 ctypes.c_int,
                                 ctypes.c_int]

def zernikeGrid(np_array, scale, m, n, radius = None, center = None):
    if (np_array.shape[0] != np_array.shape[1]):
        print("Array must be square.")
        return
    
    if radius is None:
        radius = np_array.shape[0]/2
        
    if center is None:
        #center = (np_array.shape[0]/2 - 0.5)
        center = np_array.shape[0] * 0.5

    c_np_array = numpy.ascontiguousarray(np_array, dtype=numpy.float64)

    zernike.zernike_grid(c_np_array,
                         np_array.shape[0],
                         center,
                         radius,
                         scale,
                         m,
                         n)
    
    return c_np_array


if (__name__ == "__main__"):
    print(zernike.zernike(1, 13, 0.12345, 0.0))
