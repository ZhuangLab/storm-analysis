#!/usr/bin/python
#
# Does Voroni cluster analysis.
#
# Hazen 09/16
#

import os

import storm_analysis.dbscan.cluster_stats as clusterStats
import storm_analysis.dbscan.cluster_size as clusterSize

import storm_analysis.voronoi.voronoi as voronoi

def voronoiAnalysis(bin_file, density_factor, output_directory, min_size = 30):

    # save a record of the clustering parameters.    
    bin_dir = os.path.dirname(bin_file)
    if (len(bin_dir) == 0):
        bin_dir = "."
    
    with open(bin_dir + "/voroni.txt", "w") as fp:
        fp.write("density factor = " + str(density_factor) + "\n")
        fp.write("min_size = " + str(min_size) + "\n")


    cl_bin_file = output_directory + os.path.basename(bin_file)[:-8] + "srt_list.bin"

    # find clusters
    if True:
        voronoi.voronoi(bin_file, cl_bin_file, density_factor, min_size)

    # cluster stats
    if True:
        clusterStats.clusterStats(cl_bin_file, min_size - 1)

    # cluster size
    if True:
        clusterSize.clusterSize(cl_bin_file, cl_bin_file[:-8] + "size_list.bin")


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Voronoi based clustering following Levet, Nature Methods, 2015')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations input file. This is a binary file in Insight3 format.")
    parser.add_argument('--density', dest='density', type=float, required=True,
                        help = "The density multiplier to be in a cluster. The median polygon size is multiplied by this value to give a threshold for polygon size to be in a cluster.")
    parser.add_argument('--dir', dest='output_dir', type=str, required=True,
                        help = "The directory to save the clustering results files in.")
    parser.add_argument('--min_size', dest='min_size', type=int, required=False, default=30,
                        help = "The minimum cluster size to include when calculating cluster statistics.")

    args = parser.parse_args()

    voronoiAnalysis(args.mlist, args.density, args.output_dir, args.min_size)


