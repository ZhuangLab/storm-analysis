
Getting started:

$ python ./dbscan_analysis.py --help

$ python ./cluster_images.py --help


DBSCAN analysis will create 4 output files:

dbscan.txt - A record of what the Voronoi parameters were.

X_clusters_list.bin - A molecule list file with the molecules cluster id stored in the "lk" field.

X_clusters_size_list.bin - Same as above, but in addition the cluster size is stored in the "a" field and the cluster id field is also stored in the "fr" field.

X_clusters_stats.txt - A text file containing some statistics for each of the clusters.


Note:

(1) Molecules that are not assigned to a cluster will have an "lk" value of -1. Cluster numbering starts at 2.

(2) In the default configuration clustering is done in 2D and the molecules category is ignored. To change this you will have to edit the dbscan_analysis.py or find_clusters.py

