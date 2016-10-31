
Python files:

batch_analysis.py - Run spline analysis on multiple movies in the same
   directory.

cubic_fit_c.py - The Python wrapper of the C spline fitting library.

cubic_spline_c.py - The Python wrapper of the C spline library. Given
   the splines coefficients this can be used to calculate the spline
   and it's derivatives.

find_peaks_fista.py - Peak finding using a compressed sensing approach
   and the FISTA solver.

find_peaks_std.py - Peak finding using a similar approach to 3D-DAOSTORM,
   start with the brightest local maxima and the proceed to dimmer
   maxima.

measure_psf.py - Measure the PSF of the microscope given a z stack movie
   and the location of the emitters in each frame.
   
measure_psf_beads.py - This is similar to measure_psf except that it uses
   a text for input rather than a molecule list file (.bin file).

offset_to_Z.py - Converts a .off file (from storm-control) to a z_offset
   file that can be used for measuring the PSF.

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

multi_fit_core.c - C library for spline fitting. This is compiled
   together with cubic_fit to create the C library that does the fitting.


Sample run (in the spliner/sample_data directory):

1. Take a z stack of beads -> beads_zcal.tif

2. Create a text file containing the locations of a few (isolated)
   beads -> beads_locs.txt

3. Create a text file containing the z value of each frame in the z
   stack movie -> beads_zoffset.text

4. Use measure_psf_beads.py to measure the PSF:

$ python
>>> import storm_analysis.spliner.measure_psf_beads as measure_psf_beads
>>> measure_psf_beads.measurePSFBeads("beads_zcal.tif", "beads_zoffset.txt", "beads_locs.txt", "beads_psf.psf")
('Processing frame:', 0)
('Processing frame:', 50)
('Processing frame:', 100)
('Processing frame:', 150)
('Processing frame:', 200)
('Processing frame:', 250)

   Note:
    The default measurement range is z +- 600nm in 50nm steps.
   
   This will create two files:
    (1) beads_psf.psf - The PSF file.
    (2) psf_bead.dax (and .inf) - The measured PSF as a z-stack.

5. Use psf_to_spline.py to convert the measured PSF into a spline that can be
   used by spliner/cspline for analyzing STORM movies:

>>> import storm_analysis.spliner.psf_to_spline as psf_to_spline
>>> psf_to_spline.psfToSpline("beads_psf.psf", "beads_psf.spline", 12)
Generating 3D spline.
Generating XY splines.
Generating fitting spline.
Calculating spline coefficients.

   Note:
    "12" is the size of the spline in pixels in x and y.

   This will create two files:
    (1) beads_psf.spline - A file containing the spline coefficients.
    (2) spline.dax - The spline as a z-stack. This is 2x upsampled so the its
        final size is 24 x 24 x 24.

6. Use spline_analysis.py to analyze a STORM movie. In this example this is
   just another bead movie.

>>> import storm_analysis.spliner.spline_analysis as spline_analysis
>>> spline_analysis.analyze("beads_test.tif", "beads_slist.bin", "spline_fit.xml")
Peak finding
(' Removing negative values in frame', 0)
('Frame:', 0, 23, 23)
(' Removing negative values in frame', 1)
('Frame:', 1, 23, 46)

... etc ...

(' Removing negative values in frame', 48)
('Frame:', 48, 23, 1127)
(' Removing negative values in frame', 49)
('Frame:', 49, 23, 1150)

('Added', 1150)

Tracking
Molecules: 1150 (beads_slist.bin)
Descriptor: 2
Processing molecule 0 in frame 0 (tracker)
Finished processing
Found 1150 tracks
Analysis complete


   Note:
    Most of the parameters in spline_fit.xml file are pretty similar to those in
    3D-DAOSTORM. You will need to set the "spline" parameter to the spline that
    you want to use for fitting. You will likely also need to adjust the
    "threshold" parameter depending on your data.


Optional:

If you want to refine the spline model of the PSF you can use the above spline
to bootstrap.

1. Use spline_analysis.py to analyze the calibration movie:

>>> spline_analysis.analyze("beads_zcal.tif", "beads_zcal_slist.bin", "spline_fit.xml")
Peak finding
(' Removing negative values in frame', 0)
('Frame:', 0, 16, 16)
(' Removing negative values in frame', 1)
('Frame:', 1, 16, 32)

... etc ...

(' Removing negative values in frame', 298)
('Frame:', 298, 16, 4521)
(' Removing negative values in frame', 299)
('Frame:', 299, 16, 4537)

('Added', 4537)

Tracking
Molecules: 4537 (beads_zcal_slist.bin)
Descriptor: 2
Processing molecule 0 in frame 0 (tracker)
Finished processing
Found 4537 tracks
Analysis complete

2. Use measure_psf.py to measure the PSF:

>>> import storm_analysis.spliner.measure_psf as measure_psf
>>> measure_psf.measurePSF("beads_zcal.tif", "beads_zoffset.txt", "beads_zcal_slist.bin", "beads_psf_2.psf")
('Version:', 'M425')
('Frames:', 1)
('Status:', 6)
('Molecules:', 4537)

Measuring 3D PSF
Using z offset file.
(0, 'peaks in', 16, ', peaks out', 12)
(1, 'peaks in', 16, ', peaks out', 12)

... etc ...

(298, 'peaks in', 16, ', peaks out', 12)
(299, 'peaks in', 16, ', peaks out', 12)
(0, 0.0)
(1, 0.0)

... etc ...

(29, 0.0)
(30, 0.0)

   Note:
    If you use a filename that does not exist ("foo" for example) for the z_offset file
    then fit z value in the .bin file we be used instead.

3. Use psf_to_spline.py to convert the measured PSF into a spline that can be
   used by spliner/cspline for analyzing STORM movies:

>>> psf_to_spline.psfToSpline("beads_psf_2.psf", "beads_psf_2.spline", 12)
Generating 3D spline.
Generating XY splines.
Generating fitting spline.
Calculating spline coefficients.
