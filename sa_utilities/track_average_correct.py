#!/usr/bin/python
#
# Does tracking, averaging and drift correction on a molecule list file.
#
# Hazen 10/12
#

import numpy
import os
import random
import subprocess
import sys

import sa_library.parameters as params
import sa_library.readinsight3 as readinsight3
import sa_library.writeinsight3 as writeinsight3


src_dir = os.path.dirname(__file__)

# Get z range
def getZRange(parameters):
    min_z = -0.5
    max_z = 0.5
    if hasattr(parameters, "min_z"):
        min_z = parameters.min_z
    if hasattr(parameters, "max_z"):
        max_z = parameters.max_z
    return [min_z, max_z]

# Averages all the molecules in a track into a single molecule.
def averaging(mol_list_filename, ave_list_filename):
    proc_params = [src_dir + "/avemlist",
                   mol_list_filename,
                   ave_list_filename]
    subprocess.call(proc_params)

# Performs drift correction.
def driftCorrection(mlist_file, parameters):
    drift_path = src_dir + "/"
    drift_name = mlist_file[:-4] + "_drift.txt"
    proc_params = ["python",
                   drift_path + "xyz-drift-correction.py",
                   mlist_file,
                   drift_name,
                   str(parameters.frame_step),
                   str(parameters.d_scale)]
    subprocess.call(proc_params)
    proc_params = [drift_path + "apply-drift-correction",
                   mlist_file,
                   drift_name]
    subprocess.call(proc_params)

# Does the frame-to-frame tracking.
def tracking(mol_list_filename, parameters):
    [min_z, max_z] = getZRange(parameters)
    proc_params = [src_dir + "/tracker",
                   mol_list_filename,
                   parameters.descriptor,
                   str(parameters.radius),
                   str(1000.0*min_z),
                   str(1000.0*max_z)]
    subprocess.call(proc_params)


# Peform analysis if called from the command line

if __name__ == "__main__":

    # setup
    if(len(sys.argv)==4):
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        parameters = params.Parameters(sys.argv[3])
    else:
        print "usage: track_average_correct.py <input_list.bin> <output_list.bin> <params>"
        exit()

    # tracking
    print "Tracking"
    tracking(input_file, parameters)

    # averaging
    print "Averaging"
    alist_file = input_file
    if(parameters.radius > 0.0):
        alist_file = output_file
        averaging(input_file, alist_file)
    print ""

    # drift correction
    print "Drift Correction"
    if hasattr(parameters, "drift_correction"):
        if parameters.drift_correction:
            driftCorrection(alist_file, parameters)
    print ""


#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
