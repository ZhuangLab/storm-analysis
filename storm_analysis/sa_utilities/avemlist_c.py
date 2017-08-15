#!/usr/bin/env python
"""
Python interface to the C avemlist library.

Hazen 10/16
"""

import ctypes
import os
from xml.etree import ElementTree

from storm_analysis import asciiString

import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.readinsight3 as readinsight3

c_avemlist = loadclib.loadCLibrary("storm_analysis.sa_utilities", "avemlist")

c_avemlist.avemlist.argtypes = [ctypes.c_int,
                                ctypes.c_void_p]

def avemlist(input_filename, output_filename):

    # Load input file meta data (if any).
    meta_data = readinsight3.loadI3Metadata(input_filename)
    
    argc = 3
    argv = (ctypes.c_char_p * argc)()
    argv[:] = [asciiString(elt) for elt in ["avemlist",
                                            input_filename,
                                            output_filename]]
    c_avemlist.avemlist(argc, argv)

    # Add the same meta data to the output file.
    if meta_data is not None:
        with open(output_filename, 'ab') as fp:
            fp.write(ElementTree.tostring(meta_data, 'ISO-8859-1'))
        

if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = 'Create averaged localization file from a (tracked) localization file.')

    parser.add_argument('--input_bin', dest='input_bin', type=str, required=True,
                        help = "The name of the (tracked) localization file.")
    parser.add_argument('--output_bin', dest='output_bin', type=str, required=True,
                        help = "The name of the averaged localization file.")

    args = parser.parse_args()
    
    avemlist(args.input_bin, args.output_bin)

