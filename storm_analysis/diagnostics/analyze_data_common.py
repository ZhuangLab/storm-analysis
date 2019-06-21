#!/usr/bin/env python
"""
That which is common to all of the XYZ.analyze_data diagnostic modules.

Hazen 06/21
"""
import glob
import os
import time


def analyzeData(analyzer, xml_file):
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
        analyzer.analyze(a_dir + "/test.tif", mlist, xml_file)
        stop_time = time.time()

        # Save timing results.
        print("Analysis completed in {0:.2f} seconds".format(stop_time - start_time))

        with open(a_dir + "/timing.txt", "w") as fp:
            fp.write(str(stop_time - start_time) + "\n")
