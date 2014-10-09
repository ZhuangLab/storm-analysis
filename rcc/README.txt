
The RCC code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/
Python - matplotlib - http://matplotlib.org/
(64 bit windows installers for these libraries are available here http://www.lfd.uci.edu/~gohlke/pythonlibs/)

You will need to add the project root directory (i.e. the directory one
up from the location of this file) to your Python path.

Files:
 rcc-drift-correction.py - The program for calculating drift during STORM
   movies using the RCC approach.

Usage:
 > python rcc-drift-correction.py movie_mlist.bin movie_drift.txt 500 2

The file movie_mlist.bin is the input molecule list file. The results are
saved in the movie_drift.txt file which specifies and the calculated drift
in x,y and z for every frame. 500 is the number of frames to bin into a
single STORM image for the purposes of calculating correlations between
frames. 2 is the up-scaling factor to use for the STORM image. For example,
if the original data is 256x256 then the STORM image will be 512x512.

Also, if you add an additional parameter on the end then the z drift
calculation will not be performed:

 > python rcc-drift-correction.py movie_mlist.bin movie_drift.txt 500 2 1
