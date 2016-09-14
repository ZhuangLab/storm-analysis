
The DBSCAN code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/
(64 bit windows installers for these libraries are available here http://www.lfd.uci.edu/~gohlke/pythonlibs/)

randomcolor-py - randomColor for Python - https://github.com/kevinwuhoo/randomcolor-py


A sample run:

python ./dbscan_analysis.py molecule_list.bin channel_number


This will create 3 output files:

molecule_clusters_list.bin - A molecule list file with the molecules cluster stored in the "lk" field.

molecule_clusters_size_list.bin - Same as above, but in addition the cluster size is stored in the "a" field and the cluster id field is also stored in the "fr" field.

molecule_clusters_stats.txt - A text file containing some statistics for each of the clusters.


Note:

(1) Molecules that are not assigned to a cluster will have an "lk" value of -1. Cluster numbering starts at 2.

(2) To change the DBSCAN parameters (epsilon, minimum counts and minimum cluster size) you need to edit the dbscan_analysis.py file.

(3) In the default configuration clustering is done in 2D and the molecules category is ignored.

