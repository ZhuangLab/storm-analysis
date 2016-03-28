#!/usr/bin/env python
#
# Run the 3D-DAOSTORM and sCMOS tests. You must run this in
# the test directory in order for it to work.
#
# Hazen 03/2016
#

import os
import subprocess

dao_exe = "../3d_daostorm/mufit_analysis.py"
scmos_exe = "../sCMOS/scmos_analysis.py"

test_args = [[dao_exe, "test.dax", "test_3d_2d_fixed.bin", "test_3d_2d_fixed.xml"],
             [dao_exe, "test_low_snr.dax", "test_3d_2d_fixed_low_snr.bin", "test_3d_2d_fixed_low_snr.xml"],
             [dao_exe, "test.dax", "test_3d_2d.bin", "test_3d_2d.xml"],
             [dao_exe, "test.dax", "test_3d_3d.bin", "test_3d_3d.xml"],
             [dao_exe, "test.dax", "test_3d_Z.bin", "test_3d_Z.xml"],
             [scmos_exe, "test.dax", "test_sc_2d_fixed.bin", "test_sc_2d_fixed.xml"],
             [scmos_exe, "test.dax", "test_sc_2d.bin", "test_sc_2d.xml"],
             [scmos_exe, "test.dax", "test_sc_3d.bin", "test_sc_3d.xml"],
             [scmos_exe, "test.dax", "test_sc_Z.bin", "test_sc_Z.xml"]]


# Remove any old bin files.
for arg in test_args:
    try:
        os.remove(arg[2])
    except OSError:
        pass
        

# Run analysis
for arg in test_args:
    proc_params = ["python"] + arg
    subprocess.call(proc_params)
    print ""
    print "-----"
    print ""

