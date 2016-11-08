
Gettings started:

$ python ./voronoi_analysis.py --help


Voronoi analysis will create 4 output files:

voronoi.txt - A record of what the Voronoi parameters were.

X_srt_list.bin - A molecule list file with the molecules cluster id stored in the "lk" field.

X_srt_size_list.bin - Same as above, but in addition the cluster size is stored in the "a" field and the cluster id field is also stored in the "fr" field.

X_srt_stats.txt - A text file containing some statistics for each of the clusters.


Note:

(1) Molecules that are not assigned to a cluster will have an "lk" value of -1. Cluster numbering starts at 2.

(2) The density factor is relative to the median Voronoi polygon area.

(3) Clustering is done in 2D and the molecules category is ignored.

