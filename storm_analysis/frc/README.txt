
The FRC code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/
Python - matplotlib - http://matplotlib.org/
(64 bit windows installers for these libraries are available here http://www.lfd.uci.edu/~gohlke/pythonlibs/)


In order to use the FRC code you will first need to compile the frc.c C helper 
library. Compilation instructions are included in the header of this file. A 
pre-compiled library is provided for 64-bit Windows.

Also, you will need to add the project root directory (i.e. the directory one
up from the location of this file) to your Python path.

Files:
 frc.c - The C helper library.
 frc_c.py - A Python wrapper of frc.c
 frc_calc2d.py - The program for calculating the 2D FRC.

Usage:
 > python frc_calc2d.py movie_list.bin movie_frc.txt

Results are plotted using matplotlib and also saved in movie_frc.txt.

Note that this is the un-corrected FRC, so repeated localizations of the same
molecule could artificially increase the apparent resolution.
