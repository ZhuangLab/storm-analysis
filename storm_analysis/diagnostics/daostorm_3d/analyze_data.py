#!/usr/bin/env python
"""
Analyze test data using 3D-DAOSTORM.

Hazen 09/17
"""
import glob
import os
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
    
        # Remove stale results, if any.
        if os.path.exists(hdf5):
            os.remove(hdf5)

        # Run analysis.
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
    
