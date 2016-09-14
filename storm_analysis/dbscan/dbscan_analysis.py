#!/usr/bin/python
#
# Does DBSCAN cluster analysis.
#
# Hazen 01/12
#

import os
import subprocess
import sys

if (len(sys.argv) != 3):
    print("usage <bin_file> <channel>")
    exit()

src_dir = os.path.dirname(__file__)
if not (src_dir == ""):
    src_dir += "/"

# setup
bin_file = sys.argv[1]
channel = sys.argv[2]

# defaults
eps = 40
mc = 10
min_size = 50

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
if 1:
    subprocess.call(['python', find_clusters_exe, bin_file, str(eps), str(mc)])
cl_bin_file = bin_file[:-8] + "clusters_list.bin"

# cluster stats
if 1:
    subprocess.call(['python', cluster_stats_exe, cl_bin_file, str(min_size-1)])

# cluster size
if 1:
    subprocess.call(['python', cluster_size_exe, cl_bin_file, cl_bin_file[:-8] + "size_list.bin"])

