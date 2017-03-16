.. highlight:: none

Analysis Programs
=================

These are the different localization finding and fitting approaches
in this project.

3D-DAOSTORM
-----------

This approach performs maximum likelihood estimation (MLE) guassian fitting.
It can be used to analyze 2D and 3D astigmatism STORM movies.

``storm-analysis/storm_analysis/daostorm/``

Usage ::

  # Python
  >>> from storm_analysis.daostorm_3d.mufit_analysis import analyze
  >>> analyze("movie_01.dax", "movie_01_mlist.bin", "analysis_params.xml")

  # Command line
  $ python path/to/mufit_analysis.py --movie movie_01.tif --bin movie_01_mlist.bin --xml analysis_params.xml

.. note:: The `_mlist.bin` extension is an important part of the name and
	  it is best not to substitute this for something else.   
     
Ref - `Babcock et al <http://dx.doi.org/10.1186/2192-2853-1-6>`_

sCMOS
-----

This is many ways very similar to 3D-DAOSTORM, but it is designed to handle
the analysis of data from sCMOS cameras. In order for this to work well
you will need to have a calibration file containing the offset, gain
and variance for each camera pixel.

``storm-analysis/storm_analysis/sCMOS``

Usage ::

  # Python
  >>> from storm_analysis.sCMOS.scmos_analysis import analyze
  >>> analyze("movie_01.dax", "movie_01_mlist.bin", "analysis_params.xml")

  # Command line
  $ python path/to/scmos_analysis.py --movie movie_01.tif --bin movie_01_mlist.bin --xml analysis_params.xml
  
Ref - `Huang et al <http://dx.doi.org/10.1038/nmeth.2488>`_

Spliner
-------

This approach performs MLE fitting using a cubic spline approximation of
the microscope PSF. It can be used to analyze both 2D and 3D STORM movies
with arbitrary PSF shapes. In order to use it you will need to have
a fairly accurate measurement of your microscope PSF.

``storm-analysis/storm_analysis/spliner``

Measuring the PSF
~~~~~~~~~~~~~~~~~

1. Take a z stack of beads, ``beads_zcal.tif``

2. Create a text file containing the locations of a few (isolated)
   beads, ``beads_locs.txt``

3. Create a text file containing the z value of each frame in the z
   stack movie, ``beads_zoffset.text``

4. Use measure_psf_beads.py to measure the PSF. ::

     $ python
     >>> import storm_analysis.spliner.measure_psf_beads as measure_psf_beads
     >>> measure_psf_beads.measurePSFBeads("beads_zcal.tif", "beads_zoffset.txt", "beads_locs.txt", "beads_psf.psf")
   
   This will create two files:
	  
   * beads_psf.psf - The PSF file.
  
   * psf_bead.dax (and .inf) - The measured PSF as a z-stack.

   .. note:: The default measurement range is z +- 600nm in 50nm steps.
	  
	     The default AOI size is 12 pixels.
     
Converting the PSF to a spline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

Use *psf_to_spline.py* to convert the measured PSF into a spline that can be
used by spliner for analyzing STORM movies. ::

  $ python
  >>> import storm_analysis.spliner.psf_to_spline as psf_to_spline
  >>> psf_to_spline.psfToSpline("beads_psf.psf", "beads_psf.spline", 12)

.. note:: "12" is the size of the spline in pixels in x and y.

This will create two files:

* beads_psf.spline - A file containing the spline coefficients.
  
* spline.dax - The spline as a z-stack. This is 2x upsampled so (in this examples) its
  final size is 24 x 24 x 24.

Running spliner
~~~~~~~~~~~~~~~

Edit the *analysis_params.xml* file to use *beads_psf.spline* as the spline for analysis. ::
   
   # Python
   >>> import storm_analysis.spliner.spline_analysis as spline_analysis
   >>> spline_analysis.analyze("movie_01.tif", "movie_01_slist.bin", "analysis_params.xml")

   # Command line
   $ python path/to/spline_analysis.py --movie movie_01.tif --bin movie_01_slist.bin --xml analysis_params.xml

.. note:: We use `_slist.bin` as the extension to distinguish the results from those
	  from 3D-DAOSTORM / sCMOS.

Optional
~~~~~~~~

You can refine the spline model of the PSF by using the spline determined as above to bootstrap. ::

  # Run spliner on the bead file.
  >>> spline_analysis.analyze("beads_zcal.tif", "beads_zcal_slist.bin", "analysis_params.xml")

  # Re-measure the PSF.
  >>> import storm_analysis.spliner.measure_psf as measure_psf
  >>> measure_psf.measurePSF("beads_zcal.tif", "beads_zoffset.txt", "beads_zcal_slist.bin", "beads_psf_2.psf")

  # Generate the refined spline.
  >>> psf_to_spline.psfToSpline("beads_psf_2.psf", "beads_psf_2.spline", 12)

Ref - `Babcock and Zhuang <http://dx.doi.org/10.1101/083402>`_

L1H
---

This is a compressed sensing approach. It is substantially slower than
all of the above approaches and only works with 2D STORM movies. If your
localization density is very high it may be a better choice.

``storm-analysis/storm_analysis/L1H``

Usage ::
  
  # python
  >>> from storm_analysis.L1H.cs_analysis import analyze
  >>> analyze("movie_01.dax", "movie_01.xml", "movie_01.hres", "movie_01_cslist.bin")

Ref - `Babcock et al <http://dx.doi.org/10.1364/OE.21.028583>`_
