#!/usr/bin/env python
#
# Run tests of various components of this project.
#
# Hazen 03/2016
#

import os
import subprocess

dao_exe = "../daostorm_3d/mufit_analysis.py"
l1h_exe = "../L1H/cs_analysis.py"
scmos_exe = "../sCMOS/scmos_analysis.py"

test_args = [["Testing 3D-DAOSTORM"],
             [dao_exe, "test.dax", "test_3d_2d_fixed.bin", "test_3d_2d_fixed.xml"],
             [dao_exe, "test_low_snr.dax", "test_3d_2d_fixed_low_snr.bin", "test_3d_2d_fixed_low_snr.xml"],
             [dao_exe, "test.dax", "test_3d_2d.bin", "test_3d_2d.xml"],
             [dao_exe, "test.dax", "test_3d_3d.bin", "test_3d_3d.xml"],
             [dao_exe, "test.dax", "test_3d_Z.bin", "test_3d_Z.xml"],
             ["Testing L1H"],
             ["../L1H/setup_A_matrix.py", "theoritical", "test_l1h", "1.0"],
             [l1h_exe, "test_l1h.dax", "test_l1h.xml", "test_l1h.hres", "test_l1h_list.bin"],
             ["Testing sCMOS"],
             [scmos_exe, "test.dax", "test_sc_2d_fixed.bin", "test_sc_2d_fixed.xml"],
             [scmos_exe, "test.dax", "test_sc_2d.bin", "test_sc_2d.xml"],
             [scmos_exe, "test.dax", "test_sc_3d.bin", "test_sc_3d.xml"],
             [scmos_exe, "test.dax", "test_sc_Z.bin", "test_sc_Z.xml"]]


#
# Clean up from old tests.
#

# Create a list of files to remove. old results files.
to_remove = []
for arg in test_args:
    if (len(arg) != 1):

        # Remove 3D-DAOSTORM or sCMOS results.
        if (arg[0] == dao_exe) or (arg[0] == scmos_exe):
            to_remove.append(arg[2])

        # Remove L1H results.
        elif (arg[0] == l1h_exe):
            to_remove.append("test_l1h_a7_k5_i8_o8_p8_4.amat")
            to_remove.append(arg[3])
            to_remove.append(arg[4])


# Remove the files (if the exist)
for a_file in to_remove:
    try:
        os.remove(a_file)
    except OSError:
        pass


#
# Run new tests.
#
for arg in test_args:
    if (len(arg) == 1):
        print("")
        print(arg[0])
        print("")
    else:
        proc_params = ["python3"] + arg
        subprocess.call(proc_params)
        print("")
        print("-----")
        print("")

