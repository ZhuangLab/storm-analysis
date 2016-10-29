#!/usr/bin/python
#
# Does Voroni cluster analysis.
#
# Hazen 09/16
#

import os
import subprocess
import sys


def voronoiAnalysis(bin_file, density_factor, output_directory, min_size = 30):
    src_dir = os.path.dirname(__file__)
    if not (src_dir == ""):
        src_dir += "/"
    
    bin_dir = os.path.dirname(bin_file)
    if (len(bin_dir) == 0):
        bin_dir = "."
    
    with open(bin_dir + "/voroni.txt", "w") as fp:
        fp.write("density factor = " + str(density_factor) + "\n")
        fp.write("min_size = " + str(min_size) + "\n")

    # exe files
    voroni_exe = src_dir + "/voronoi.py"
    cluster_stats_exe = src_dir + "../dbscan/cluster_stats.py"
    cluster_size_exe = src_dir + "../dbscan/cluster_size.py"

    cl_bin_file = output_directory + os.path.basename(bin_file)[:-8] + "srt_list.bin"

    # find clusters
    if True:
        subprocess.call(['python', voroni_exe, bin_file, str(density_factor), str(min_size), cl_bin_file])

    # cluster stats
    if True:
        subprocess.call(['python', cluster_stats_exe, cl_bin_file, str(min_size-1)])

    # cluster size
    if True:
        subprocess.call(['python', cluster_size_exe, cl_bin_file, cl_bin_file[:-8] + "size_list.bin"])


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Voronoi based clustering following Levet, Nature Methods, 2015')

    parser.add_argument('--bin', dest='mlist', type=str, required=True)
    parser.add_argument('--density', dest='density', type=float, required=True)
    parser.add_argument('--dir', dest='output_dir', type=str, required=True)
    parser.add_argument('--min_size', dest='min_size', type=int, required=False, default=30)

    args = parser.parse_args()

    voronoiAnalysis(args.mlist, args.density, args.output_dir, args.min_size)


