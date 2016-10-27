#!/usr/bin/env python
#
# Python interface to the C tracker library. Note that this
# library uses static variables so it is not thread safe.
#
# Hazen 10/16
#

import ctypes
import os

import storm_analysis.sa_library.loadclib as loadclib

ctracker = loadclib.loadCLibrary(os.path.dirname(__file__), "tracker")

ctracker.tracker.argtypes = [ctypes.c_int,
                             ctypes.c_void_p]

def tracker(mlist_filename, descriptor, radius, zmin, zmax, save_track_id = 0):
    argc = 7
    argv = (ctypes.c_char_p * argc)()
    argv[:] = ["tracker",
               mlist_filename,
               descriptor,
               str(radius),
               str(zmin),
               str(zmax),
               str(save_track_id)]
    ctracker.tracker(argc, argv)

if (__name__ == "__main__"):
    import sys
    
    tracker(*sys.argv[1:])
