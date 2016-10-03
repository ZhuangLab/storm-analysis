#!/usr/bin/python
#
# Perform spline analysis on a dax file given parameters.
#
# Hazen 01/16
#

import sys

import storm_analysis.spliner.find_peaks_fista as find_peaks_fista
import storm_analysis.spliner.find_peaks_std as find_peaks_std
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as std_analysis

# setup
if(len(sys.argv)==3):
    parameters = params.Parameters(sys.argv[2])
    mlist_file = sys.argv[1][:-4] + "_mlist.bin"
elif(len(sys.argv)==4):
    parameters = params.Parameters(sys.argv[3])
    mlist_file = sys.argv[2]
else:
    print("usage: <movie> <bin> <parameters.xml>")
    exit()

if hasattr(parameters, "use_fista") and (parameters.use_fista != 0):
    finder = find_peaks_fista.SplinerFinderFitter(parameters)
else:
    finder = find_peaks_std.SplinerFinderFitter(parameters)
    
std_analysis.standardAnalysis(finder,
                              sys.argv[1],
                              mlist_file,
                              parameters)


