
Files:

astigmaticPSF.py - Describes an astigmatic PSF.

background.py - Simulation background generation classes.

camera.py - Simulation camera emulator classes.

check_photophysics.py - For verifying that the (truth) emitter photophysics is as expected.

dhPSF.py - Describes a double-helix PSF.

draw_gaussians.c
draw_gaussians_c.py - Draws gaussian shaped peaks on an image.

emitters_on_grid.py - Create i3.bin files with emitters on a X/Y grid.

photophysics.py - Simulation dye photophysics classes.

psf.py - Simulation emitter rendering classes.

pupil_math.py - PSF generation using a pupil function approach.

simbase.py - Base class for the various simulation classes.

simulate.py - Generates simulated data using classes such camera, psf, etc. that specify the simulation behavior.

simulate_2d.py - Very simple 2D simulator.

simulate_3d.py - Very simple 3D simulator.

zernike.c - C library for (more) rapid calculation of zernike polynomials.
zernike_c.py - Python interface to the C zernike library.


Note:

Please see the simulate.py file for an explanation of the recommended way to perform simulations, as well as a simple example. The goal is to break the simulation process into interchangeable modules for each step in the simulation process.
