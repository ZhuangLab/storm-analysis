#!/usr/bin/python
#
# Does Voroni cluster analysis.
#
# Hazen 09/16
#

import os
import subprocess
import sys

if (len(sys.argv) != 4):
    print("usage <bin_file> <density factor> <output directory>")
    exit()

src_dir = os.path.dirname(__file__)
if not (src_dir == ""):
    src_dir += "/"
    
# setup
bin_file = sys.argv[1]
output_directory = sys.argv[3]

# defaults
density_factor = float(sys.argv[2])
min_size = 50

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
if 1:
    subprocess.call(['python', voroni_exe, bin_file, str(density_factor), cl_bin_file])

# cluster stats
if 1:
    subprocess.call(['python', cluster_stats_exe, cl_bin_file, str(min_size-1)])

# cluster size
if 1:
    subprocess.call(['python', cluster_size_exe, cl_bin_file, cl_bin_file[:-8] + "size_list.bin"])

