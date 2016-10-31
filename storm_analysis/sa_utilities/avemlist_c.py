#!/usr/bin/env python
#
# Python interface to the C avemlist library.
#
# Hazen 10/16
#

import ctypes
import os

from storm_analysis import asciiString
import storm_analysis.sa_library.loadclib as loadclib

c_avemlist = loadclib.loadCLibrary("storm_analysis.sa_utilities", "_avemlist")

c_avemlist.avemlist.argtypes = [ctypes.c_int,
                                ctypes.c_void_p]

def avemlist(input_filename, output_filename):
    argc = 3
    argv = (ctypes.c_char_p * argc)()
    argv[:] = [asciiString(elt) for elt in ["avemlist",
                                            input_filename,
                                            output_filename]]
    c_avemlist.avemlist(argc, argv)

if (__name__ == "__main__"):
    import sys
    
    avemlist(*sys.argv[1:])
