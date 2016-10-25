
The cubic spline code has the same dependencies as 3D-DAOSTORM and sCMOS.


Python files:

batch_analysis - Run spline analysis on multiple movies in the same
   directory.

cubic_fit_c - The Python wrapper of the C spline fitting library.

cubic_spline_c - The Python wrapper of the C spline library. Given
   the splines coefficients this can be used to calculate the spline
   and it's derivatives.

find_peaks_fista - Peak finding using a compressed sensing approach
   and the FISTA solver.

find_peaks_std - Peak finding using a similar approach to 3D-DAOSTORM,
   start with the brightest local maxima and the proceed to dimmer
   maxima.

measure_psf - Measure the PSF of the microscope given a z stack movie
   and the location of the emitters in each frame.
   
measure_psf_beads - This is similar to measure_psf except that it uses
   a text for input rather than a molecule list file (.bin file).

offset_to_Z - Converts a .off file (from storm-control) to a z_offset
   file that can be used for measuring the PSF.

psf_to_spline - Generate a spline that describes the measured PSF.

spline1D - 1D cubic spline in Python.

spline2D - 2D cubic spline in Python.

spline3D - 3D cubic spline in Python.

spline_analysis - Run spline analysis on a single STORM movie.

spline_to_psf - Generate a PSF at a particular z value given a spline.


C files:

cubic_fit - C library for spline fitting.

cubic_spline - C library for generating the PSF and it's deratives
   from a cubic spline.

multi_fit_core - C library for spline fitting. This is compiled
   together with cubic_fit to create the C library that does the fitting.


Sample Run:

These commands should be executed in the spliner/sample_data directory:

1. Take a z stack of beads -> beads_zcal.tif

2. Create a text file containing the locations of a few (isolated)
   beads -> beads_locs.txt

3. Create a text file containing the z value of each frame in the z
   stack movie -> beads_zoffset.text

4. Use measure_psf_beads.py to measure the PSF:

   > python ../measure_psf_beads.py beads_zcal.tif beads_zoffset.txt beads_locs.txt beads_psf.psf

   Note:
   The default measurement range is z +- 600nm in 50nm steps.
   
   This will create two files:
   (1) beads_psf.psf - The PSF file.
   (2) psf_bead.dax (and .inf) - The measured PSF as a z-stack.

5. Use psf_to_spline.py to convert the measured PSF into a spline that can be
   used by spliner/cspline for analyzing STORM movies:

   > python ../psf_to_spline.py beads_psf.psf beads_psf.spline 12

   Note:
   "12" is the size of the spline in pixels in x and y.

   This will create two files:
   (1) beads_psf.spline - A file containing the spline coefficients.
   (2) spline.dax - The spline as a z-stack. This is 2x upsampled so the its
       final size is 24 x 24 x 24.

6. Use spline_analysis.py to analyze a STORM movie. In this example this is
   just another bead movie.

   > python ../spline_analysis.py beads_test.tif beads_slist.bin spline_fit.xml

   Note:
   Most of the parameters in spline_fit.xml file are pretty similar to those in
   3D-DAOSTORM. You will need to set the "spline" parameter to the spline that
   you want to use for fitting. You will likely also need to adjust the
   "threshold" parameter depending on your data.


Optional:

If you want to refine the spline model of the PSF you can use the above spline
to bootstrap.

1. Use spline_analysis.py to analyze the calibration movie:

   > python ../spline_analysis.py beads_zcal.tif beads_zcal_slist.bin spline_fit.xml

2. Use measure_psf.py to measure the PSF:

   > python ../measure_psf.py beads_zcal.tif beads_zoffset.txt beads_zcal_slist.bin beads_psf_2.psf 1

   Notes:
   (1) "1" - This tells measure_psf.py that we want a 3D PSF (as opposed to 2D).
   (2) If you use a filename that does not exist ("foo" for example) for the z_offset file
       then fit z value in the .bin file we be used instead.

3. Use psf_to_spline.py to convert the measured PSF into a spline that can be
   used by spliner/cspline for analyzing STORM movies:

   > python ../psf_to_spline.py beads_psf_2.psf beads_psf_2.spline 12   
