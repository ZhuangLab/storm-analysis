#!/usr/bin/env python
"""
Analyze SLURM test data using Multiplane. We don't actually use
SLURM, instead we just analyze each movie sequentially.

Hazen 09/18
"""
import glob
import os
import time

import storm_analysis.multi_plane.multi_plane as mp

def analyzeData():
    jobs = sorted(glob.glob("slurm_test/job*.xml"))

    for job in jobs:
        print()
        print("Analyzing job:", job)
        print()
        
        h5_name = os.path.splitext(job)[0] + ".hdf5"

        # Remove stale results, if any.
        if os.path.exists(h5_name):
            os.remove(h5_name)

        # Run analysis.
        start_time = time.time()
        mp.analyze("slurm_test/test", h5_name, job)
        stop_time = time.time()

        # Print timing results.
        print("Analysis completed in {0:.2f} seconds".format(stop_time - start_time))


if (__name__ == "__main__"):
    analyzeData()
    
