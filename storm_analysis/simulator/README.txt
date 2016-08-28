
The simulator code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/
(64 bit windows installers for these libraries are available here http://www.lfd.uci.edu/~gohlke/pythonlibs/)


Files:

astigmaticPSF.py - Describes an astigmatic PSF.

compile_linux.sh - Run this to compile the C library on linux.

dhPSF.py - Describes a double-helix PSF.

draw_gaussians.c
drawgaussians.py - Draws gaussian shaped peaks on an image.

simulate_2d.py - Creates simulated 2D STORM images.

simulate_3d.py - Creates simulated 3D STORM images.


Note:

These are very simple simulations. The molecules are either on a grid or randomly
distributed. The noise is Poisson and the background is flat.
