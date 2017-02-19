#!/usr/bin/env python
"""
Does DBSCAN cluster analysis.

Hazen 01/12
"""

import os

import storm_analysis.dbscan.find_clusters as findClusters
import storm_analysis.dbscan.cluster_stats as clusterStats
import storm_analysis.dbscan.cluster_size as clusterSize


def dbscanAnalysis(bin_file, channel, eps = 40, mc = 10, min_size = 50):

    # save a record of the clustering parameters.
    bin_dir = os.path.dirname(bin_file)
    if (len(bin_dir) == 0):
        bin_dir = "."
    
    with open(bin_dir + "/dbscan.txt", "w") as fp:
        fp.write("eps = " + str(eps) + "\n")
        fp.write("mc = " + str(mc) + "\n")
        fp.write("min_size = " + str(min_size) + "\n")

    cl_bin_file = bin_file[:-8] + "clusters_list.bin"

    # find clusters
    if True:
        findClusters.findClusters(bin_file, cl_bin_file, eps, mc)

    # cluster stats
    if True:
        clusterStats.clusterStats(cl_bin_file, min_size - 1)

    # cluster size
    if True:
        clusterSize.clusterSize(cl_bin_file, cl_bin_file[:-8] + "size_list.bin")


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'DBSCAN clustering following Ester, KDD-96, 1996')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations input file. This is a binary file in Insight3 format.")
    parser.add_argument('--channel', dest='channel', type=int, required=True,
                        help = "Which channel (or category) to use for clustering.")
    parser.add_argument('--eps', dest='epsilon', type=float, required=False, default=40,
                        help = "The DBSCAN epsilon parameters in nanometers. The default is 40nm.")
    parser.add_argument('--mc', dest='mc', type=int, required=False, default=10,
                        help = "The DBSCAN mc parameter. The default is 10.")
    parser.add_argument('--min_size', dest='min_size', type=int, required=False, default=50,
                        help = "The minimum cluster size to include when calculating cluster statistics. The default is 50.")

    args = parser.parse_args()

    dbscanAnalysis(args.mlist, args.channel, args.epsilon, args.mc, args.min_size)
