#!/usr/bin/env python
"""
Simple Python interface to pf_math.c

Hazen 01/19
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.loadclib as loadclib

pf_math = loadclib.loadCLibrary("pf_math")

pf_math.pfZernikeGrid.argtypes = [ndpointer(dtype=numpy.float64),
                                  ctypes.c_int,
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_int,
                                  ctypes.c_int]


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
        # This center matches that of pupil_math.Geometry.
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
