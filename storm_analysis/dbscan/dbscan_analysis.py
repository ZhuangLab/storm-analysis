#!/usr/bin/python
#
# Does DBSCAN cluster analysis.
#
# Hazen 01/12
#

import os
import subprocess
import sys


def dbscanAnalysis(bin_file, channel, eps = 40, mc = 10, min_size = 50):
    src_dir = os.path.dirname(__file__)
    if not (src_dir == ""):
        src_dir += "/"

    bin_dir = os.path.dirname(bin_file)
    if (len(bin_dir) == 0):
        bin_dir = "."
    
    with open(bin_dir + "/dbscan.txt", "w") as fp:
        fp.write("eps = " + str(eps) + "\n")
        fp.write("mc = " + str(mc) + "\n")
        fp.write("min_size = " + str(min_size) + "\n")

    # exe files
    find_clusters_exe = src_dir + "find_clusters.py"
    cluster_stats_exe = src_dir + "cluster_stats.py"
    cluster_size_exe = src_dir + "cluster_size.py"

    # find clusters
    if True:
        subprocess.call(['python', find_clusters_exe, bin_file, str(eps), str(mc)])
    cl_bin_file = bin_file[:-8] + "clusters_list.bin"

    # cluster stats
    if True:
        subprocess.call(['python', cluster_stats_exe, cl_bin_file, str(min_size-1)])

    # cluster size
    if True:
        subprocess.call(['python', cluster_size_exe, cl_bin_file, cl_bin_file[:-8] + "size_list.bin"])


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'DBSCAN clustering following Ester, KDD-96, 1996')

    parser.add_argument('--bin', dest='mlist', type=str, required=True)
    parser.add_argument('--channel', dest='channel', type=int, required=True)
    parser.add_argument('--eps', dest='epsilon', type=float, required=False, default=40)
    parser.add_argument('--mc', dest='mc', type=int, required=False, default=10)
    parser.add_argument('--min_size', dest='min_size', type=int, required=False, default=50)

    args = parser.parse_args()

    dbscanAnalysis(args.mlist, args.channel, args.epsilon, args.mc, args.min_size)
