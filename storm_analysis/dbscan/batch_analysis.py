#!/usr/bin/env python
"""
Batch clustering analysis.

Hazen 01/12
"""

import glob
import os

import storm_analysis.sa_library.batch_run as batchRun


def batchAnalysis(input_directory, eps = 40, mc = 10, min_size = 50):
    """
    Batch DBSCAN clustering.
    """
    # We use the executable and batchRun as this is an easy
    # way to run all the clustering in parallel.
    src_dir = os.path.dirname(__file__)
    if not (src_dir == ""):
        src_dir += "/"
    
    clusters_exe = src_dir + "dbscan_analysis.py"

    # Find appropriate bin files.
    bin_files = glob.glob(input_directory + "*.hdf5")

    # Generate command line.
    cmd_lines = []
    for filename in bin_files:

        # skip clustering related bin files
        if ("clusters" in filename) or ("srt" in filename):
            continue

        print("Found:", filename)

        cmd_lines.append(['python', clusters_exe,
                          "--bin", filename,
                          "--eps", str(eps),
                          "--mc", str(mc),
                          "--min_size", str(min_size)])

    batchRun.batchRun(cmd_lines)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Batch DBSCAN clustering.')

    parser.add_argument('--dir', dest='directory', type=str, required=True,
                        help = "The name of the directory containing the Insight3 format files to analyze.")
    parser.add_argument('--eps', dest='epsilon', type=float, required=False, default=40,
                        help = "The DBSCAN epsilon parameters in nanometers. The default is 40nm.")
    parser.add_argument('--mc', dest='mc', type=int, required=False, default=10,
                        help = "The DBSCAN mc parameter. The default is 10.")
    parser.add_argument('--min_size', dest='min_size', type=int, required=False, default=50,
                        help = "The minimum cluster size to include when calculating cluster statistics. The default is 50.")

    args = parser.parse_args()

    batchAnalysis(args.directory, args.channel, args.epsilon, args.mc, args.min_size)

