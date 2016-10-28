#!/usr/bin/env python
#
# Python interface to the C apply_drift_correction library. Note that this
# library uses static variables so it is not thread safe.
#
# Hazen 10/16
#

import ctypes
import os

import storm_analysis.sa_library.loadclib as loadclib

adc = loadclib.loadCLibrary("storm_analysis.sa_utilities", "_apply-drift-correction")

adc.applyDriftCorrection.argtypes = [ctypes.c_int,
                                     ctypes.c_void_p]

def applyDriftCorrection(mlist_filename, drift_filename):
    argc = 3
    argv = (ctypes.c_char_p * argc)()
    argv[:] = ["apply-drift-correction",
               mlist_filename,
               drift_filename]
    adc.applyDriftCorrection(argc, argv)

if (__name__ == "__main__"):
    import sys
    
    applyDriftCorrection(*sys.argv[1:])
