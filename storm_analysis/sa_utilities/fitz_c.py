#!/usr/bin/env python
#
# Python interface to the C fitz library. Note that this
# library uses static variables so it is not thread safe.
#
# Hazen 10/16
#

import ctypes
import os

import storm_analysis.sa_library.loadclib as loadclib

c_fitz = loadclib.loadCLibrary("storm_analysis.sa_utilities", "_fitz")

c_fitz.fitz.argtypes = [ctypes.c_int,
                        ctypes.c_void_p]

def fitz(i3_filename, cut_off, wx_params, wy_params):
    argc = 17
    argv = (ctypes.c_char_p * argc)()
    argv[0] = "fitz"
    argv[1] = i3_filename
    argv[2] = str(cut_off)
    for i in range(7):
        argv[i+3] = str(wx_params[i])
        argv[i+10] = str(wy_params[i])
    c_fitz.fitz(argc, argv)

if (__name__ == "__main__"):
    import sys

    if(len(sys.argv) != 17):
        print("fitz_c.py requires 16 arguments")
        exit()
    
    fitz(sys.argv[1], sys.argv[2], sys.argv[3:10], sys.argv[10:17])
