#!/usr/bin/python
#
# Does tracking, averaging and drift correction on a molecule list file.
#
# Tracking is performed "in place" on the input_list.bin file.
# Averaging is creates (or overwrites) the output_list.bin file based
#    on the input_list.bin file.
# Drift correction is performed in place one or both file depending
#    on whether or not averaging was done.
#
# Hazen 10/13
#

import os
import sys

import sa_library.parameters as params
import sa_utilities.std_analysis as std_analysis

# Setup
if(len(sys.argv)==4):
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    parameters = params.Parameters(sys.argv[3])
else:
    print "usage: track_average_correct.py <input_list.bin> <output_list.bin> <params>"
    exit()

# Tracking
print "Tracking"
std_analysis.tracking(input_file, parameters)

# Averaging
print "Averaging"
did_averaging = False
if(parameters.radius > 0.0):
    did_averaging = True
    std_analysis.averaging(input_file, output_file)
print ""

# Drift correction
print "Drift Correction"
if hasattr(parameters, "drift_correction"):
    if parameters.drift_correction:
        if did_averaging:
            std_analysis.driftCorrection([input_file, output_file], parameters)
        else:
            std_analysis.driftCorrection([input_file], parameters)
print ""

#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
