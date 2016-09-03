#!/usr/bin/env python
#
# Run tests of various components of this project. Depending
# on the speed of your computer this could take a few minutes
# to run.
#
# Hazen 08/2016
#

import os
import subprocess

dao_exe = "../daostorm_3d/mufit_analysis.py"
scmos_exe = "../sCMOS/scmos_analysis.py"
spliner_exe = "../spliner/spline_analysis.py"

test_args = [["Testing 3D-DAOSTORM"],
             [dao_exe, "test.dax", "test_3d_2d_fixed.bin", "test_3d_2d_fixed.xml"],
             [dao_exe, "test_low_snr.dax", "test_3d_2d_fixed_low_snr.bin", "test_3d_2d_fixed_low_snr.xml"],
             [dao_exe, "test.dax", "test_3d_2d.bin", "test_3d_2d.xml"],
             [dao_exe, "test.dax", "test_3d_3d.bin", "test_3d_3d.xml"],
             [dao_exe, "test.dax", "test_3d_Z.bin", "test_3d_Z.xml"],
             ["Testing L1H"],
             ["../L1H/setup_A_matrix.py", "theoritical", "test_l1h", "1.0", "0"],
             ["../L1H/cs_analysis.py", "test_l1h.dax", "test_l1h.xml", "test_l1h.hres", "test_l1h_list.bin"],
             ["Testing Rolling Ball Background Subtraction"],
             ["../rolling_ball_bgr/rolling_ball.py", "test_bg_sub.dax", "test_bg_sub_rb.dax", "10", "1"],
             ["Testing sCMOS"],
             [scmos_exe, "test.dax", "test_sc_2d_fixed.bin", "test_sc_2d_fixed.xml"],
             [scmos_exe, "test.dax", "test_sc_2d.bin", "test_sc_2d.xml"],
             [scmos_exe, "test.dax", "test_sc_3d.bin", "test_sc_3d.xml"],
             [scmos_exe, "test.dax", "test_sc_Z.bin", "test_sc_Z.xml"],
             ["Testing Spliner"],
             ["../spliner/measure_psf.py", "test_spliner.dax", "none", "test_spliner_olist.bin", "test_spliner_psf.psf", "1"],
             ["../spliner/psf_to_spline.py", "test_spliner_psf.psf", "test_spliner_psf.spline", "10"],
             [spliner_exe, "test_spliner.dax", "test_spliner_slist.bin", "test_spliner_dh.xml"],
             [spliner_exe, "test_spliner.dax", "test_spliner_flist.bin", "test_spliner_dh_fista.xml"],
             ["Testing Wavelet Background Subtraction"],
             ["../wavelet_bgr/wavelet_bgr.py", "test_bg_sub.dax", "test_bg_sub_wbgr.dax", "db4", "2", "2", "10"]]


#
# Clean up from old tests.
#

# Create a list of files to remove.
to_remove = []

# 3D-DAOSTORM and sCMOS results.
for arg in test_args:
    if (len(arg) != 1):
        if (arg[0] == dao_exe) or (arg[0] == scmos_exe):
            to_remove.append(arg[2])

# L1H results.
to_remove.append("test_l1h_a7_k5_i8_o8_p8_4.amat")
to_remove.append("test_l1h.hres")
to_remove.append("test_l1h_list.bin")

# Rolling ball background subtraction results
to_remove.append("test_bg_sub_rb.dax")
to_remove.append("test_bg_sub_rb.inf")

# Spliner results.
to_remove.append("test_spliner_psf.psf")
to_remove.append("test_spliner_psf.spline")
to_remove.append("test_spliner_slist.bin")
to_remove.append("test_spliner_alist.bin")
to_remove.append("test_spliner_flist.bin")

# Wavelet background subtraction results
to_remove.append("test_bg_sub_wbgr.dax")
to_remove.append("test_bg_sub_wbgr.inf")


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

