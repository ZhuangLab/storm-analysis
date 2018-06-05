
Python files:

batch_analysis.py - Run spline analysis on multiple movies in the same
   directory.

cramer_rao.py - Calculate the Cramer-Rao bound given a 3D spline.

cubic_fit_c.py - The Python wrapper of the C spline fitting library.

cubic_spline_c.py - The Python wrapper of the C spline library. Given
   the splines coefficients this can be used to calculate the spline
   and it's derivatives.

find_peaks_fista.py - Peak finding using a compressed sensing approach
   and the FISTA solver.

find_peaks_std.py - Peak finding using a similar approach to 3D-DAOSTORM,
   start with the brightest local maxima and the proceed to dimmer
   maxima.

hdf5_to_beads.py - Convert an HDF5 file to a beads file for PSF
   measurement.

measure_psf.py - Measure the PSF of the microscope given a z stack movie
   and the location of the emitters in each frame.
   
measure_psf_beads.py - This is similar to measure_psf except that it uses
   a text for input rather than a molecule list file (.bin file).

measure_psf_utils.py - Utility functions used for PSF measurement.

offset_to_Z.py - Converts a .off file (from storm-control) to a z_offset
   file that can be used for measuring the PSF.

print_psf.py - Print out some information about a PSF file.

psf_to_spline.py - Generate a spline that describes the measured PSF.

spline1D.py - 1D cubic spline in Python.

spline2D.py - 2D cubic spline in Python.

spline3D.py - 3D cubic spline in Python.

spline_analysis.py - Run spline analysis on a single STORM movie.

spline_to_psf.py - Generate a PSF at a particular z value given a spline.


C files:

cubic_fit.c - C library for spline fitting.

cubic_spline.c - C library for generating the PSF and it's deratives
   from a cubic spline.


