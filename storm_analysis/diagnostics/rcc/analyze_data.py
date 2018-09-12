#!/usr/bin/env python
"""
Analyze drift correction data using RCC.

Hazen 01/18
"""
import glob
import os
import shutil
import time

import storm_analysis.rcc.rcc_drift_correction as rccDriftCorrection

import storm_analysis.diagnostics.rcc.settings as settings


def analyzeData():
    dirs = sorted(glob.glob("test*"))
    total_time = 0.0
    for a_dir in dirs:
        print()
        print("Analyzing:", a_dir)
        print()
    
        hdf5 = os.path.join(a_dir, "test.hdf5")

        # Do tracking and drift correction using the ground truth positions.
        if True:
            shutil.copyfile(os.path.join(a_dir, "test_ref.hdf5"), hdf5)

        # Run analysis.
        start_time = time.time()
        rccDriftCorrection.rccDriftCorrection(hdf5,
                                              os.path.join(a_dir, "test_drift.txt"),
                                              settings.step,
                                              settings.scale,
                                              settings.z_min,
                                              settings.z_max,
                                              settings.correct_z,
                                              make_plots = False)
        stop_time = time.time()

        # Save timing results.
        total_time += stop_time - start_time
        print("Analysis completed in {0:.2f} seconds.".format(stop_time - start_time))

        with open(os.path.join(a_dir, "timing.txt"), "w") as fp:
            fp.write(str(stop_time - start_time) + "\n")

    print()
    print("{0:d} directories analyzed in {1:.2f} seconds.".format(len(dirs), total_time))


if (__name__ == "__main__"):
    analyzeData()
    
