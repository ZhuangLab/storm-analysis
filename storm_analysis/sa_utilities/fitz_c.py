#!/usr/bin/env python
"""
Python interface to the C fitz library.

Hazen 10/16
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

from storm_analysis import asciiString
import storm_analysis.sa_library.loadclib as loadclib

c_fitz = loadclib.loadCLibrary("storm_analysis.sa_utilities", "fitz")

c_fitz.fitz.argtypes = [ctypes.c_char_p,
                        ndpointer(dtype=numpy.float64),
                        ndpointer(dtype=numpy.float64),
                        ctypes.c_double,
                        ctypes.c_double,
                        ctypes.c_double,
                        ctypes.c_double]

def fitz(i3_filename, cutoff, wx_params, wy_params, z_min, z_max, z_step = 1.0):
    """
    This expects all z related parameters to be in nanometers.
    """
    c_fitz.fitz(asciiString(i3_filename),
                numpy.ascontiguousarray(wx_params),
                numpy.ascontiguousarray(wy_params),
                cutoff,
                z_min,
                z_max,
                z_step)
                

if (__name__ == "__main__"):
    
    import argparse

    import storm_analysis.sa_library.parameters as params

    parser = argparse.ArgumentParser(description = 'Z fitting given Wx, Wy calibration curves.')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations file. This is a binary file in Insight3 format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()

    parameters = params.ParametersDAO().initFromFile(args.settings)

    [wx_params, wy_params] = parameters.getWidthParams()
    [min_z, max_z] = parameters.getZRange()
        
    fitz(args.mlist,
         parameters.getAttr("cutoff"),
         wx_params,
         wy_params,
         min_z * 1000.0,
         max_z * 1000.0,
         parameters.getAttr("z_step", 1.0))
