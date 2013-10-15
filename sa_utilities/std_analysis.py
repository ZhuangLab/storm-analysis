#!/usr/bin/python
#
# Performs "standard" analysis on a dax file given parameters.
#
# Hazen 10/13
#

import os
import subprocess
import sys

import sa_library.parameters as params

src_dir = os.path.dirname(__file__)
if not (src_dir == ""):
    src_dir += "/"

# Averages all the molecules in a track into a single molecule.
def averaging(mol_list_filename, ave_list_filename):
    proc_params = [src_dir + "avemlist",
                   mol_list_filename,
                   ave_list_filename]
    subprocess.call(proc_params)

# Performs drift correction.
def driftCorrection(list_files, parameters):
    drift_name = list_files[0][:-9] + "drift.txt"

    # Check if we have been asked not to do z drift correction.
    # The default is to do the correction.
    z_correct = True
    if hasattr(parameters, "z_correction"):
        if not parameters.z_correction:
            z_correct = False

    if z_correct:
        proc_params = ["python",
                       src_dir + "xyz-drift-correction.py",
                       list_files[0],
                       drift_name,
                       str(parameters.frame_step),
                       str(parameters.d_scale)]
    else:
        proc_params = ["python",
                       src_dir + "xyz-drift-correction.py",
                       list_files[0],
                       drift_name,
                       str(parameters.frame_step),
                       str(parameters.d_scale),
                       str(1)]
    subprocess.call(proc_params)

    if (os.path.exists(drift_name)):
        for list_file in list_files:
            proc_params = [src_dir + "apply-drift-correction",
                           list_file,
                           drift_name]
            subprocess.call(proc_params)

# Perform standard analysis.
def standardAnalysis(peakFindingFn, data_file, mlist_file, parameters):

    # peak finding
    print "Peak finding"
    if(not peakFindingFn(data_file, mlist_file, parameters)):
        print ""
        
        # tracking
        print "Tracking"
        tracking(mlist_file, parameters)

        # averaging
        alist_file = None
        if(parameters.radius > 0.0):
            alist_file = mlist_file[:-9] + "alist.bin"
            averaging(mlist_file, alist_file)
            print ""

        # z fitting
        if parameters.do_zfit:
            print "Fitting Z"
            if alist_file:
                zFitting(alist_file, parameters)
            zFitting(mlist_file, parameters)
            print ""

        # drift correction
        if hasattr(parameters, "drift_correction"):
            if parameters.drift_correction:
                print "Drift Correction"
                if alist_file:
                    driftCorrection([mlist_file, alist_file], parameters)
                else:
                    driftCorrection([mlist_file], parameters)
                print ""
    print "Analysis complete"

# Does the frame-to-frame tracking.
def tracking(mol_list_filename, parameters):
    [min_z, max_z] = params.getZRange(parameters)
    proc_params = [src_dir + "tracker",
                   mol_list_filename,
                   parameters.descriptor,
                   str(parameters.radius),
                   str(1000.0*min_z),
                   str(1000.0*max_z)]
    subprocess.call(proc_params)

# Does z fitting.
def zFitting(mol_list_filename, parameters):
    if(parameters.orientation == "inverted"):
        wx_str = map(str, params.getWidthParams(parameters, "y"))
        wy_str = map(str, params.getWidthParams(parameters, "x"))
    else:
        wx_str = map(str, params.getWidthParams(parameters, "x"))
        wy_str = map(str, params.getWidthParams(parameters, "y"))
    proc_params = [src_dir + "fitz",
                   mol_list_filename, 
                   str(parameters.cutoff)] + wx_str + wy_str
    subprocess.call(proc_params)

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
