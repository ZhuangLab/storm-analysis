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
``storm_analysis/simulator/draw_gaussians_c.py`` to render the
localizations as gaussians onto a numpy array. The *getPSFs()*
method of the *GaussianPSF* class in ``storm_analysis/simulator/psf.py``
provides an example of how to do this.

Clustering
----------

DBSCAN
~~~~~~

Voronoi
~~~~~~~

Other
-----

FRC
~~~

