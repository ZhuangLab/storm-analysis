.. highlight:: none

Analysis Programs
=================

These are the different localization finding and fitting approaches in this project.

Note that all of the fitting based approaches (i.e. everything except L1H) follow the
approach described in `Tang et al <http://dx.doi.org/10.1038/srep11073>`_ for localization
identification. This means that you will need to know how to convert from ADU to
photo-electrons for your camera in order to get the best results.

All of the fitting approaches use the Levenberg-Marquadt algorithm for maximum likelihood
estimation fitting (MLE) and achieve the Cramer-Rao bound for localization accuracy.

3D-DAOSTORM
-----------

This approach performs maximum likelihood estimation (MLE) Gaussian fitting.
It can be used to analyze 2D and 3D astigmatism STORM movies.

``storm-analysis/storm_analysis/daostorm_3d/``

Usage ::

  # Python
  >>> from storm_analysis.daostorm_3d.mufit_analysis import analyze
  >>> analyze("movie_01.dax", "movie_01.hdf5", "analysis_params.xml")

  # Command line
  $ python path/to/mufit_analysis.py --movie movie_01.tif --bin movie_01.hdf5 --xml analysis_params.xml

``storm-analysis/storm_analysis/daostorm_3d/z_calibration.py`` can used to measure
the Z fitting parameters necessary for 3D astigmatism imaging.

Ref - `Babcock et al <http://dx.doi.org/10.1186/2192-2853-1-6>`_

sCMOS
-----

This is essentially identical to 3D-DAOSTORM, but it is designed to handle
the analysis of data from sCMOS cameras. In order for this to work well
you will need to have a calibration file containing the offset, gain, variance
and relative QE for each camera pixel.

``storm-analysis/storm_analysis/sCMOS``

Usage ::

  # Python
  >>> from storm_analysis.sCMOS.scmos_analysis import analyze
  >>> analyze("movie_01.dax", "movie_01.hdf5", "analysis_params.xml")

  # Command line
  $ python path/to/scmos_analysis.py --movie movie_01.tif --bin movie_01.hdf5 --xml analysis_params.xml
  
Ref - `Huang et al <http://dx.doi.org/10.1038/nmeth.2488>`_

Spliner
-------

This approach performs MLE fitting using a cubic spline approximation of
the microscope PSF. It can be used to analyze both 2D and 3D STORM movies
with arbitrary PSF shapes. In order to use it you will need to have
a fairly accurate measurement of your microscope PSF.

``storm-analysis/storm_analysis/spliner``

It accepts either EMCCD or sCMOS camera data.

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

   .. note:: The default measurement range is z +- 0.6um in 0.05um steps.
	  
	     The default value for `aoi_size` is 12 pixels. The actual size of the
	     measured AOI will be 2x this value, or 24 pixels in this example.
     
Converting the PSF to a spline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

Use *psf_to_spline.py* to convert the measured PSF into a spline that can be
used by spliner for analyzing STORM movies. ::

  $ python
  >>> import storm_analysis.spliner.psf_to_spline as psf_to_spline
  >>> psf_to_spline.psfToSpline("beads_psf.psf", "beads_psf.spline", 10)

.. note:: "10" is 1/2 the size of the spline in pixels in x and y.

This will create two files:

* beads_psf.spline - A file containing the spline coefficients.
  
* spline.dax - The spline as a z-stack. The final size will be 20 x 20 x 20
  in this example as the "size" of the spline would be more correctly termed
  1/2 the size of the spline.

Running spliner
~~~~~~~~~~~~~~~

Edit the *analysis_params.xml* file to use *beads_psf.spline* as the spline for analysis. ::
   
   # Python
   >>> import storm_analysis.spliner.spline_analysis as spline_analysis
   >>> spline_analysis.analyze("movie_01.tif", "movie_01.hdf5", "analysis_params.xml")

   # Command line
   $ python path/to/spline_analysis.py --movie movie_01.tif --bin movie_01.hdf5 --xml analysis_params.xml


Optional
~~~~~~~~

You can refine the spline model of the PSF by using the spline determined as above to bootstrap. ::

  # Run spliner on the bead file.
  >>> spline_analysis.analyze("beads_zcal.tif", "beads_zcal.hdf5", "analysis_params.xml")

  # Re-measure the PSF.
  >>> import storm_analysis.spliner.measure_psf as measure_psf
  >>> measure_psf.measurePSF("beads_zcal.tif", "beads_zoffset.txt", "beads_zcal.hdf5", "beads_psf_2.psf")

  # Generate the refined spline.
  >>> psf_to_spline.psfToSpline("beads_psf_2.psf", "beads_psf_2.spline", 10)

Ref - `Babcock and Zhuang <http://dx.doi.org/10.1101/083402>`_


Multiplane
-----------

This approach performs MLE fitting using a cubic spline approximation of the microscope PSF for
multiplane (and single plane) sCMOS data. It can be used to analyze 3D STORM movies with arbitrary
PSF shapes. In order to use it you will need to have a fairly accurate measurement of your microscope
PSF as well as transforms between the different planes.

Multiplane assumes that you have a separate movie for each channel. In what follows we will assume
that the first channel movie is called ``movie_01_ch1.tif``, the second is ``movie_01_ch2.tif`` and
etc...

If the movies are from different cameras the cameras are expected to be synchronized, i.e. they are all
exposing at the same time, and they are not all free running independently of each other. It is okay
however if they don't agree on the frame number as this can be compensated for with the
``channelX_offset`` parameter.

.. note:: Most of the scripts referenced below are in ``storm-analysis/storm_analysis/multi_plane`` folder.
	  All of them are in the ``storm-analysis`` project.
	  
``storm-analysis/storm_analysis/multi_plane``

Camera sCMOS calibration
~~~~~~~~~~~~~~~~~~~~~~~~

You will need one sCMOS calibration file per channel/plane. These are the same format as used in
the sCMOS analysis package described above.

Plane to plane mapping
~~~~~~~~~~~~~~~~~~~~~~

Multiplane analysis requires information about how to map localization XY positions in one channel
to XY positions in another channel. This can be done using the following steps:

1. Acquire a movie with reasonably bright, small and well separated beads, ``map_01_ch1.tif``,
   ``map_01_ch2.tif``, etc.. If there is a large z separation between the planes you may need
   to scan the focus during the movie.

2. Analyze one frame of each channel with sCMOS or possibly 3D-DAOSTORM to localize the beads,
   ``map_01_ch1.hdf5``, ``map_01_ch2.hdf5``, etc.. For each channel you probably
   want one of the frames that is in focus.

3. Identify the mappings between ch1 and the other channels using micrometry. ::
	  
      # Command line
      $ python path/to/micrometry/micrometry.py --locs1 map_01_ch1.hdf5 --locs2 map_01_ch2.hdf5 --results c1_c2_map.map
      $ python path/to/micrometry/micrometry.py --locs1 map_01_ch1.hdf5 --locs2 map_01_ch3.hdf5 --results c1_c3_map.map
      $ ..

   .. note:: You may need to change the ``--max_size`` parameter (in pixels) depending on how sparse your beads are.

   .. note:: You can also use the PyQt5 GUI program ``mapper.py`` for determining the channel to channel maps.
	     
4. Merge the individual mapping files using merge_maps.py. ::
	  
      # Command line
      $ python path/to/micrometry/merge_maps.py --results map.map --maps c1_c2_map.map c1_c3_map.map c1_c4_map.map ...

   .. note:: The individual mapping files must be listed in the channel order, lowest to highest.    

Measuring the PSFs
~~~~~~~~~~~~~~~~~~

1. Take a z stack of beads, ``beads_zcal_ch1.tif``, ``beads_zcal_ch2.tif``, etc..

2. Analyze one frame of the channel 1 bead movie with sCMOS or possibly 3D-DAOSTORM to localize
   the beads, ``beads_zcal_ch1.hdf5``.

3. Select good localizations to use for PSF determination for each channel. ::

     # Command line
     $ python path/to/psf_localizations.py --bin beads_zcal_ch1.hdf5 --map map.map --aoi_size 12

   .. note:: An AOI size of 12 pixels is appropriate for setups with a camera pixel size of ~100nm.

4. Create averaged z stacks for each channel. ::

     # Command line
     $ python path/to/psf_zstack.py --movie beads_zcal_ch1.tif --bin beads_zcal_ch1_c0_psf.hdf5 --zstack ch1_zstack --scmos_cal ch1_cal.npy --aoi_size 12
     $ python path/to/psf_zstack.py --movie beads_zcal_ch2.tif --bin beads_zcal_ch1_c1_psf.hdf5 --zstack ch2_zstack --scmos_cal ch2_cal.npy --aoi_size 12
     $ ..

   .. note:: (Linear) drift during the PSF calibration movie can be adjusted for using the
	     ``--driftx`` and ``--drifty`` parameters. Units are pixels per frame.
   
   .. note:: Drift can be estimated with the program ``zstack_xydrift.py``. You will need to
	     have found localizations in the first and last frame of the PSF calibration movie.

5. Create a text file containing the z offset of each frame of the PSF calibration movie. One
   possibility is to use ``spliner/offset_to_z.py``.

6. Measure the PSF for each plane. ::

     # Command line
     $ python path/to/measure_psf.py --zstack ch1_zstack.npy --zoffsets z_offsets.txt --psf_name ch1_psf.psf
     $ ..

   .. note:: You can adjust the z range of the PSF measurement using the ``z_range`` parameter.
   
   .. note:: At this point it is probably a good idea to check your PSF using a tool like ImageJ.
	  
   .. note:: If you are doing spectrally resolved STORM (`SR-STORM <http://dx.doi.org/10.1038/nmeth.3528>`_)
	     include the command line argument ``--normalize`` and skip the next step.

7. Normalize the PSFs relative to each other. ::
     
     # Command line
     $ python path/to/normalize_psfs.py --psfs ch1_psf.psf ch2_psf.psf ..

8. (Optional) Check plane z offsets using ``check_plane_offsets.py``. If the offsets are not well
   centered this can be adjusted using the ``--deltaz`` argument to ``spliner/offset_to_z.py`` and
   restarting at step 5.
     
Converting the PSFs to a splines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the same procedure as for ``Spliner``.

Use *psf_to_spline.py* to convert the measured PSF into a spline that can be
used by spliner for analyzing STORM movies. ::
  
  # Command line (if you used normalize_psfs.py).
  $ python path/to/spliner/psf_to_spline.py --psf ch1_psf_normed.psf --spline ch1_psf.spline --spline_size 10
  $ ..

  # Command line (if you did not use normalize_psfs.py).
  $ python path/to/spliner/psf_to_spline.py --psf ch1_psf.psf --spline ch1_psf.spline --spline_size 10
  $ ..

.. note:: Using ``--spline_size 10`` is appropriate for setups with a camera pixel size of ~100nm. The final
	  spline will be 20x20 pixels in X/Y.
	  
Creating the Weights File
~~~~~~~~~~~~~~~~~~~~~~~~~

Multiplane uses channel "information" weights in order to more optimally weight the contribution
from each plane in the determination of a localizations parameters. The channels are weighted
based on their Cramer-Rao bounds as a function of z.

1. Create a multiplane analysis XML file ``movie_01_analysis.xml``. A sample is available here:
   ``multi_plane/sample_data/example_analysis.xml``. Use a value of ``1`` for the
   ``independent_heights`` parameter when doing SR-STORM analysis.

2. Create the weights file. ::
	
     # Command line (all planes have the same background).
     $ python path/to/plane_weighting.py --background 20 --photons 4000 --output weights.npy --xml movie_01_analysis.xml
     
     # Command line (the background is different in each plane).
     $ python path/to/plane_weighting.py --background 20 18 15 etc.. --photons 4000 --output weights.npy --xml movie_01_analysis.xml

   .. note:: ``--background`` is photo-electrons per plane and ``--photons`` is the expected average
	     number of photo-electrons per localization summed over all the planes. If your camera
	     does not have a gain of 1.0 you will need to convert camera counts to photo-electrons.

Running Multiplane
~~~~~~~~~~~~~~~~~~

Once you have done all of the above you are finally ready to run multiplane analysis. ::

   # Command line
   $ python path/to/multi_plane.py --basename movie_01_ --bin movie_01.hdf5 --xml movie_01_analysis.xml

.. note:: The movie names that are loaded are the concatenation of ``basename`` and the values of
	  the ``channelX_ext`` parameters.

.. note:: The script ``find_offsets.py`` is useful for determining the frame difference, if any, between
	  movies from different cameras. This can be useful if the movies did not all start at the same time.

Post-analysis
~~~~~~~~~~~~~

Multiplane will generate a HDF5 file containing the localizations for all of the channels. At this point
you can do either or both of the following. Note however that these require that you ran the tracking with
a non-zero radius.

1. Calculate the first moment of the localization height as a function of channel number. ::
   
     # Command line
     $ python path/to/channel_color.py --bin movie_01.hdf5 --order 0 2 1 3

   .. note:: The order parameter is the order of the channels by increasing (or decreasing) wavelength.

   .. note:: This will add the fields 'height_moment' and 'height_total' to the tracks.
   
2. Use k-means clustering for color determination. ::

     # Command line
     $ python path/to/kmean_measure_codebook.py --bin movie_01.hdf5 --ndyes 2 --output movie_01_codebook.npy
     $ python path/to/kmean_classifier.py --codebook movie_01_codebook.npy --bin movie_01.hdf5

   .. note:: This will add the fields 'km_color' and 'km_distance' to the tracks.
      
   .. note:: Use the expected number of different dyes for the ndyes parameter.
	  
   .. note:: The default is to put the localizations in the top 20% in terms of distance from the category center
	     into the rejects category (category 9).

   .. note:: You can use a codebook from a different sample for classification.

Ref - `Babcock <http://dx.doi.org/doi:10.1038/s41598-018-19981-z>`_

Pupil Function
--------------

This approach performs MLE fitting using a pupil function to model the microscope PSF.

``storm-analysis/storm_analysis/pupilfn``

It accepts either EMCCD or sCMOS camera data.


PSF FFT
-------

This approach performs MLE fitting using the measured PSF and the Fast Fourier Transform (FFT)
to model the microscope PSF.

``storm-analysis/storm_analysis/psf_fft``

It accepts either EMCCD or sCMOS camera data.

Like ``Pupil Function`` it was written primarily to test our claim that (cubic) splines are the most efficient way
to represent an arbitrary microscope PSF.

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
