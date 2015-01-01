
The Wavelet background removal code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/
(64 bit windows installers for these libraries are available here http://www.lfd.uci.edu/~gohlke/pythonlibs/)
Python - pywt library - http://pywavelets.readthedocs.org/en/latest/index.html

This code can be used as a standalone background removal program as well as
a module that you can import.

A sample run in standalone mode:
python ./wavelet_bgr.py /path/to/movie wavelet_type wavelet_level iterations threshold

wavelet_type - See the pywt documentation, typically something like "db4".
wavelet_level - How many levels of wavelet decomposition to perform. The larger the
   number the less response to local changes in the background.
iterations - The number of iterations of background estimation and foreground 
   replacement to perform (see the Galloway paper).
threshold - This should probably be something like 1x to 2x the estimated noise in the 
   background.

This will create a file called subtracted.dax which contains the original movie 
minus the estimated background with a 100 offset.
