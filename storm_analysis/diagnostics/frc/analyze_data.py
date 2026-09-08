#!/usr/bin/env python
"""
Analyze FRC data.

Hazen 01/18
"""
import glob
import os

import storm_analysis.frc.frc_calc2d as frcCalc2d


def analyzeData():
    dirs = sorted(glob.glob("test*"))
    total_time = 0.0
    for a_dir in dirs:
        print()
        print("Analyzing:", a_dir)
        print()
    
        hdf5 = os.path.join(a_dir, "test.hdf5")
        frc_text = os.path.join(a_dir, "frc.txt")

        # Run FRC analysis.
        # show_plot defaults to True, and this loops over every test
        # directory, so without this it stops on a window each time.
        frcCalc2d.frcCalc2d(hdf5, frc_text, show_plot = False)

    print()


if (__name__ == "__main__"):
    analyzeData()
    
