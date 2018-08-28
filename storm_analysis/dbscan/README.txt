
The DBSCAN algorithm is described in this reference:

Martin Ester, Hans-Peter Kriegel, Jörg Sander and Xiaowei Xu,
"A density-based algorithm for discovering clusters in large spatial
databases with noise", Proceedings of the Second International
Conference on Knowledge Discovery and Data Mining, 1996.


Python programs:

batch_analysis.py - Run multiple DBSCANs at once.

cluster_images.py - Python functions for making pictures of clusters.

clusters_sa_h5py.py - For reading and writing storm-analysis format
                      HDF5 files with clustering information.

dbscan_analysis - Python functions for DBSCAN analysis.

dbscan_c - Python wrapper of the DBSCAN C library.


C programs:

dbscan.c - DBSCAN C library.

kdtree.c - A KDTree implementation that is used for efficient
           3D nearest neighbor searches.
	   
