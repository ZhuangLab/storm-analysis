#!/usr/bin/python
#
# Does DBSCAN cluster analysis.
#
# Hazen 01/12
#

import subprocess
import sys

if (len(sys.argv) != 4):
    print("usage <bin_file> <output directory> <channel>")
    exit()
    
# setup
bin_file = sys.argv[1]
output_directory = sys.argv[2]
channel = sys.argv[3]

# defaults
eps = 40
mc = 10
min_size = 50

# exe files
find_clusters_exe = sys.path[0] + "/find_clusters.py"
cluster_stats_exe = sys.path[0] + "/cluster_stats.py"
cluster_size_exe = sys.path[0] + "/cluster_size.py"

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

