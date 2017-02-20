#!/usr/bin/env python
"""
Python interface to the C tracker library. Note that this
library uses static variables so it is not thread safe.

Hazen 10/16
"""

import ctypes
import os

from storm_analysis import asciiString
import storm_analysis.sa_library.loadclib as loadclib

c_tracker = loadclib.loadCLibrary("storm_analysis.sa_utilities", "tracker")

c_tracker.tracker.argtypes = [ctypes.c_int,
                              ctypes.c_void_p]

def tracker(mlist_filename, descriptor, radius, zmin, zmax, save_track_id = 0):
    argc = 7
    argv = (ctypes.c_char_p * argc)()
    argv[:] = [asciiString(elt) for elt in ["tracker",
                                            mlist_filename,
                                            descriptor,
                                            radius,
                                            zmin,
                                            zmax,
                                            save_track_id]]
    c_tracker.tracker(argc, argv)

if (__name__ == "__main__"):
    import sys
    
    tracker(*sys.argv[1:])
