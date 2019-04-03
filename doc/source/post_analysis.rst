.. highlight:: none
	    
Post Analysis
=============

There are a few utilities for processing and rendering the results of
the analysis included in this project.

.. note:: The `Matlab STORM <https://github.com/ZhuangLab/matlab-storm>`_
	  project also includes a viewer.

Visualization & Rendering
-------------------------

Visualizer
~~~~~~~~~~

``storm_analysis/visualizer/visualizer.py``

You can use this program to evaluate the quality of the analysis and
to decide if the analysis parameters are reasonable. It provides a
simple GUI viewer in which you can load a STORM movie and a XX.hdf5
or a XX_list.bin file. It will show an overlay of the localizations that
were found in the current frame.

Rendering
~~~~~~~~~

Image can be created from .hdf5 localization files using
``storm_analysis/sa_utilities/hdf5_to_image.py``. This can be used
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
works on HDF5 format localization files.

Getting started: ::

  $ python ./dbscan_analysis.py --help

  $ python ./cluster_images.py --help

Voronoi
~~~~~~~

``storm_analysis/voronoi/voronoi_analysis.py``

This is an implementation of the SR-Tesseler clustering algorithm
(`Levet et al <http://dx.doi.org/10.1038/nmeth.3579>`_) that
works on HDF5 format localization files.

Gettings started: ::

  $ python ./voronoi_analysis.py --help

Note:

* This is a pure Python implementation and may be slow when dealing with
  more than 1M localizations.

* The density factor is relative to the median Voronoi polygon area.

* Clustering is always done in 2D and the molecules category is always ignored.

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
