.. highlight:: none
	    
Post Analysis
=============

There are a few utilities for processing and rendering the results of
the analysis included in this project.

We typically use a program called **Insight3** for viewing the final STORM
image. This program is only available by request from the Zhuang lab.
The `Matlab STORM <https://github.com/ZhuangLab/matlab-storm>`_ project
also includes a viewer.

Visualization & Rendering
-------------------------

Visualizer
~~~~~~~~~~

``storm_analysis/visualizer/visualizer.py``

You can use this program to evaluate the quality of the analysis and
to decide if the analysis parameters are reasonable. It provides a
simple GUI viewer in which you can load a STORM movie and an
XX_list.bin file. It will show an overlay of the localizations that
were found in the current frame.

There are two options for XX_list.bin files to load, but the only
difference between them is how the properties of the selected
localization are displayed. Different analysis programs use
the same property name for different purposes, as explained
in the XX_list.bin section of the **Output file** page.

Rendering
~~~~~~~~~

There are basically two options here.

The first and likely simplest is to use
``storm_analysis/sa_library/i3togrid.py`` to create a 2D histogram
of your data which you can save with a module like
`tifffile <https://pypi.python.org/pypi/tifffile>`_ or
`Pillow <https://pypi.python.org/pypi/Pillow/>`_

The second and slightly more complicated way is to use
``storm_analysis/sa_utilities/bin_to_image.py`` to render the
localizations as either a histogram or gaussians. This can be used
from the command line to create 2D .tif images (32 bit float format).
In addition it provides a number of functions that you can use to
create 2D or 3D images (returned as numpy arrays).

Clustering
----------

DBSCAN
~~~~~~

``storm_analysis/dbscan/dbscan_analysis.py``

This is an implementation of the
`DBSCAN <https://en.wikipedia.org/wiki/DBSCAN>`_ clustering algorithm that
works on XX_list.bin files.

Getting started: ::

  $ python ./dbscan_analysis.py --help

  $ python ./cluster_images.py --help


DBSCAN analysis will create 4 output files:

1. dbscan.txt - A record of what the DBSCAN parameters were.

2. X_clusters_list.bin - A molecule list file with the molecules cluster
   id stored in the "lk" field.

3. X_clusters_size_list.bin - Same as above, but in addition the cluster
   size is stored in the "a" field and the cluster id field is also stored
   in the "fr" field.

4. X_clusters_stats.txt - A text file containing some statistics for each
   of the clusters.

Note:

* Molecules that are not assigned to a cluster will have an "lk" value of
  -1. Cluster numbering starts at 2.

* In the default configuration clustering is done in 2D and the molecules
  category is ignored. To change this you will have to edit dbscan_analysis.py
  or find_clusters.py

Voronoi
~~~~~~~

``storm_analysis/voronoi/voronoi_analysis.py``

This is an implementation of the SR-Tesseler clustering algorithm
(`Levet et al <http://dx.doi.org/10.1038/nmeth.3579>`_) that
works on XX_list.bin files, 

Gettings started: ::

  $ python ./voronoi_analysis.py --help

Voronoi analysis will create 4 output files:

1. voronoi.txt - A record of what the Voronoi parameters were.

2. X_srt_list.bin - A molecule list file with the molecules cluster id stored
   in the "lk" field.

3. X_srt_size_list.bin - Same as above, but in addition the cluster size is
   stored in the "a" field and the cluster id field is also stored in the "fr" field.

4. X_srt_stats.txt - A text file containing some statistics for each of the clusters.

Note:

* This is a pure Python implementation and may be slow when dealing with
  more than 1M localizations.

* Molecules that are not assigned to a cluster will have an "lk" value of -1.
  Cluster numbering starts at 2.

* The density factor is relative to the median Voronoi polygon area.

* Clustering is done in 2D and the molecules category is ignored.

Other
-----

FRC
~~~

``storm_analysis/frc/frc_calc2d.py``

This is an implementation of the FRC method
(`Nieuwenhuizen et al <http://dx.doi.org/10.1038/nmeth.2448>`_)
for estimating the resolution of a STORM image.

Getting started: ::

  $python ./frc_calc2d.py --help

Note:
  
* Results are plotted using matplotlib and also saved in the output text file.

* This is the un-corrected FRC, so repeated localizations of the same
  molecule could artificially increase the apparent resolution.
