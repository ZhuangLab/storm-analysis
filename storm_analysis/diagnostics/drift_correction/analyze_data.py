#!/usr/bin/env python
"""
Analyze drift correction data using the standard image based
drift correction approach.

Hazen 01/18
"""
import glob
import os
import shutil
import time

import storm_analysis.daostorm_3d.mufit_analysis as dao3d


def analyzeData():
    dirs = sorted(glob.glob("test*"))
    total_time = 0.0
    for a_dir in dirs:
        print()
        print("Analyzing:", a_dir)
        print()
    
        hdf5 = a_dir + "/test.hdf5"

        # Do tracking and drift correction using the ground truth positions.
        if True:
            shutil.copyfile(a_dir + "/test_ref.hdf5", hdf5)

        # Run analysis.
        #
        # This will just do the drift correction as it will see that the
        # movie has been completely analyzed.
        #
        start_time = time.time()
        dao3d.analyze(a_dir + "/test.dax", hdf5, "dao.xml")
        stop_time = time.time()

        # Save timing results.
        total_time += stop_time - start_time
        print("Analysis completed in {0:.2f} seconds.".format(stop_time - start_time))

        with open(a_dir + "/timing.txt", "w") as fp:
            fp.write(str(stop_time - start_time) + "\n")

    print()
    print("{0:d} directories analyzed in {1:.2f} seconds.".format(len(dirs), total_time))


if (__name__ == "__main__"):
    analyzeData()
    
