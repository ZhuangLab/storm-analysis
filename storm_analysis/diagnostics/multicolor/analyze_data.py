#!/usr/bin/env python
"""
Analyze test data using Multiplane

Hazen 01/18
"""
import glob
import inspect
import os
import subprocess
import time

import storm_analysis
import storm_analysis.multi_plane.multi_plane as mp

mp_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/multi_plane/"

dirs = sorted(glob.glob("test*"))

for a_dir in dirs:
    print()
    print("Analyzing:", a_dir)
    print()
    
    mlist = a_dir + "/test.hdf5"

    # Remove stale results, if any.
    if os.path.exists(mlist):
        os.remove(mlist)

    # Run analysis.
    start_time = time.time()
    mp.analyze(a_dir + "/test", mlist, "multicolor.xml")

    # Categorize results using k-means clustering.
    #
    
    # Measure codebook.
    print("Measuring k-means codebook.")
    subprocess.call(["python", mp_path + "kmeans_measure_codebook.py",
                     "--bin", a_dir + "/test.hdf5",
                     "--ndyes", "4",
                     "--output", a_dir + "/codebook.npy"])
        
    # K-means cluster.
    print("k-means clustering.")
    subprocess.call(["python", mp_path + "kmeans_classifier.py",
                     "--bin", a_dir + "/test.hdf5",
                     "--codebook", a_dir + "/codebook.npy"])
        
    stop_time = time.time()

    # Save timing results.
    print("Analysis completed in {0:.2f} seconds".format(stop_time - start_time))

    with open(a_dir + "/timing.txt", "w") as fp:
        fp.write(str(stop_time - start_time) + "\n")
