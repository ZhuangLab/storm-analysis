
The DBSCAN code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/
(64 bit windows installers for these libraries are available here http://www.lfd.uci.edu/~gohlke/pythonlibs/)


A sample run:

python ./voronoi_analysis.py molecule_list.bin 1.25 ./


This will create 3 output files:

molecule_srt_list.bin - A molecule list file with the molecules cluster stored in the "lk" field.

molecule_srt_size_list.bin - Same as above, but in addition the cluster size is stored in the "a" field and the cluster id field is also stored in the "fr" field.

molecule_srt_stats.txt - A text file containing some statistics for each of the clusters.


Note:

(1) Molecules that are not assigned to a cluster will have an "lk" value of -1. Cluster numbering starts at 2.

(2) The density factor is relative to the median Voronoi polygon area.

(3) Clustering is done in 2D and the molecules category is ignored.

