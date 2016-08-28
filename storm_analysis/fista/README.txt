
The FISTA code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/
(64 bit windows installers for these libraries are available here http://www.lfd.uci.edu/~gohlke/pythonlibs/)

FFTW - http://www.fftw.org/

Files:

compile_linux.sh - Run this to compile the C library on linux.

fista_3d.py - A pure Python implementation of FISTA for deconvolving images in 3D (and 2D).

fista_decon_utilities.c
fista_decon_utilities_c.py - These are some utility functions for processing the FISTA deconvolved images.

fista_fft.c
fista_fft_c.py - The C version of fista_3d.py, this about 4x faster.
