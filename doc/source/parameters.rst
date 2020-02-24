.. highlight:: none
	       
Analysis Parameters
===================

These **.xml** files describe how to do the analysis. The parsing
of these files is handled by ``storm_analysis.sa_library.parameters.py``.
Example setting files can be found in ``storm_analysis/test/data``.

General
-------

These parameters are common to all of the analysis programs.

If you want to limit the analysis to a sub-section of the frame you have
two choices. (1) Specify ``aoi_radius``, ``x_center`` and ``y_center`` for a
circular analysis AOI. (2) Specify one or more of ``x_start``, ``x_stop``,
``y_start`` and ``y_stop``.

*Analysis parameters.*

* **aoi_radius** - Radius in pixels for a circular analysis AOI. ``x_center`` and ``y_center`` are also required.
  
* **max_frame** - The frame to stop analysis on, ``-1`` = analyze to the end of the film.

* **max_z** - Maximum z value for z fitting, specified in um.
  
* **min_z** - Minimum z value for z fitting, specified in um.
    
* **pixel_size** - CCD pixel size (in nm).

* **start_frame** - The frame to start analysis on, ``-1`` = start at the beginning of the film.

* **static_background_estimate** - If this is set, and set to a number greater than 0,
  then the analysis will estimate the background by using the average over this number of
  frames.

  If this is not set, or set to 0, the background is estimated separately for each frame.

* **x_center** - X center for a circular analysis AOI. ``y_center`` and ``aoi_radius`` are also required.
  
* **x_start** - X start of the analysis AOI, leave unset to start at the edge of the image.

* **x_stop** - X end of the analysis AOI, leave unset to end at the edge of the image.

* **y_center** - Y center for a circular analysis AOI. ``x_center`` and ``aoi_radius`` are also required.

* **y_start** - Y start of the analysis AOI, leave unsetto start at the edge of the image.

* **y_stop** - Y end of the analysis AOI, leave unset to end at the edge of the image.

*Tracking parameters.*

* **descriptor** - Frame descriptor string. ::
    
    0 - activation frame
    1 - non-specific frame
    2 - channel1 frame
    3 - channel2 frame
    4 - etc..

    A two-color multi-activator STORM example:
    "02110311"

* **radius** - Radius for matching peaks from frame to frame. Localizations that are closer
  than this value (in pixels) in adjacent frames (ignoring activation frames) are assumed
  to come from the same emitter and are averaged together to create a (hopefully) 
  more accurately localized emitter. If this is zero then no matching will be done.

*Drift correction parameters.*

* **d_scale** - This is the "scale" at which to render the sub-STORM images for drift
  correction. Drift correction works by creating STORM images from frame_step sized groups 
  of frames. These are rendered scaled by the d_scale parameter. For example, if
  your data is 256x256 pixels then the drift-correction will create 512x512 sub-STORM 
  images (for d_scale = 2) and then attempt to correlate these images to each other
  to calculate the drift. Using a larger d_scale value creates higher resolution 
  sub-STORM images, but they are also sparser so you might not see any improvement
  in the drift correction.
  
  2 is usually a good choice.

* **drift_correction** - Do drift correction, ``0`` = No.

* **frame_step** - Number of frames in each (drift correction) sub-STORM image, something
  like 500 or 1000 is a good choice here.

* **z_correction** - Do z drift correction, ``0`` = No. The default value for this
  parameter is ``1``, so you have to set it to zero for 2D data otherwise you will
  get an error message during drift correction.

*File conversion parameters.*

* **convert_to** - Specify what, if any, formats to convert the output HDF5 file into
  upon completion of the analysis. Options are .bin and .txt. Use a comma separated list
  if you want both. i.e. ".bin, .txt".


Fitting Based Analysis
----------------------

These parameters are common to fitting based analysis approaches, 3D-DAOSTORM, sCMOS, Spliner and Multiplane.

* **background_sigma** - This is the sigma of a 2D gaussian to convolve the data in order to estimate
  the background in pixels. Something like 8 usually works well.

* **fftw_estimate** - FFTW should estimate the best FFT plan instead of measuring which is best. This can
  help speed the analysis of short movies that are very large in XY. ``0`` = (default) FFTW will
  measure the best FFT plan. ``1`` = FFTW will estimate the best FFT plan.

* **find_max_radius** - To be a peak it must be the maximum value within this radius (in pixels).

* **fit_error_model** - Specify which fitting error model to use.

  The default model is 'MLE' (Maximum Likelihood Estimation). Other options include  
  'ALS' (Anscombe Least Squares), 'LS' (Least Squares), 'DWLS' (Data Weighted Least 
  Squares) and 'FWLS' (Fit Weighted Least Squares).
                           
  Note that most fitters only support the 'MLE' fitting error model.

* **no_fitting** - If this is non-zero then we won't do any fitting iterations. This is useful for
  testing the finder, as well as how accurately we're initializing the peak parameter values.

* **iterations** - Maximum number of iterations for new peak finding. Usually you'd use
  something like 20 to try and separate neighboring peaks. However if you are analyzing
  beads or some other bright and sparse target 1 may be a better choice as it will suppress
  the spurious splitting of a single peak into multiple peaks.

* **peak_locations** - This is used when you already know where your want fitting to
  happen, as for example in a bead calibration movie and you just want to use the
  approximate locations as inputs for fitting.

  peak_locations is a text file with the peak x, y, height and background
  values as white spaced columns (x and y positions are in pixels as
  determined using **visualizer**). ::
  
    1.0 2.0 1000.0 100.0
    10.0 5.0 2000.0 200.0
    ...
  
* **sigma** - This is the estimated sigma of the PSF in pixels, it serves several
  purposes.

  (1) It is used in most of the analysis approaches as a measure of the
      peak to peak distance at which peak fits do not substantially
      effect each other.

  (2) In most of the analysis, if two peaks are closer than this distance
      then the dimmer one will be discarded.
  
  (3) In **3D-DAOSTORM** and **sCMOS** analysis it is also used as the initial guess
      for the peak sigma.

  (4) In **3D-DAOSTORM** and **sCMOS** the peak widths are constrained to be between
      0.5x and 5x the value of sigma.

  (5) In **Multiplane DAO** the peak widths are constrained to be between
      0.5x and 3x the value of sigma.
  
* **threshold** - Ideally this is in units of sigma, as in a "x sigma event". For example
  at 3 sigma you'd expect about 0.003 false positives per pixel. Incorrect background
  estimation can however complicate things. You probably want to use a value greater than
  6.0 for most analysis. Also if your label is quite bright and you are not modelling
  your peaks that well (incorrect PSF, PSF is too small) you may need to set this higher
  to avoid getting apparently spurious low intensity peaks.


3D-DAOSTORM and sCMOS
---------------------

These parameters are common to 3D-DAOSTORM and sCMOS analysis.

* **cutoff** - Z fit cutoff (used when z is calculated later from the fit width in x and y.

* **do_zfit** - Do z fitting (or not), only relevant for "3d" fitting (see "model" parameter).

* **foreground_sigma** - This is the sigma of a 2D gaussian to convolve the image with prior to peak
  indentification. When your data has a low SNR this can help for peak finding. For optimal sensitivity
  it should be the same as the expected sigma for your peaks. If you set it to zero (or comment it out)
  then this will not be performed, which can make the analysis (very slightly) faster.  

* **model** - Model is one of "2dfixed", "2d", "3d", or "Z". ::

    2dfixed - fixed sigma 2d gaussian fitting.
    2d - variable sigma 2d gaussian fitting.
    3d - x, y sigma are independently variable, z
         will be fit after peak fitting.
    Z - x, y sigma depend on z, z is fit as part
         of peak fitting.

* **roi_size** - The fitting ROI size to use in pixels. The total number of pixels is
  roi_size * roi_size. If this is not specified then it will be calculated based from
  the sigma value and the fitting model. Basically by increasing/decreasing this you
  are trading off accuracy versus speed. A value that is 6x the largest sigma you expect
  to fit is a good compromise.

* **sigma_range** - A two element array that specifies the minimum and maximum sigma values to
  allow when fitting for the peak width. If this is not specified the default is
  [0.5 * ``sigma``, 5.0 * ``sigma``]. Only relevant for the "2d" and "3d" fitting models.

* **wx vs z parameters** - These are used for determining the localization Z position
  based on its in width in x and y (astigmatism imaging). See
  `Huang et al <http://dx.doi.org/10.1126/science.1153529>`_ for
  a more detailed explanation. Units are either nanometers or dimensionless.
            
  * wx_wo
  * wx_c
  * wx_d
  * wxA
  * wxB
  * wxC
  * wxD
     
* **wy vs z parameters** - Same as above.
       
  * wy_wo
  * wy_c
  * wy_d
  * wyA
  * wyB
  * wyC
  * wyD

* **z_value** - The starting z value for fitting. If this is not specified it defaults to 0.0.
  Units are microns.

* **z_step** - The z step size for finding the optimal z value when using the 3d model. If
  this is not specified it defaults to 1 nanometer. Units are microns.

3D-DAOSTORM
-----------

* **camera_gain** - Conversion factor to go from camera ADU to photo-electrons. Units are ADU/e-,
  so the camera ADU values will be divided by this number to convert to photo-electrons (e-).

* **camera_offset** - This what the camera reads with the shutter closed.

sCMOS
-----

* **camera_calibration** - This file contains the sCMOS calibration data for the region of
  the camera that the movie comes from. It consists of 4 numpy arrays, [offset, variance, gain,
  relative QE], each of which is the same size as a frame of the movie that is to be analyzed.
  This can be generated for a camera using camera_calibration.py and (if it needs
  to be resliced), reslice_calibration.py.

Spliner
-------

* **spline** - This is the spline file to use for fitting. Based on the spline the analysis
  will decide whether to do 2D or 3D spline fitting, 2D if the spline is 2D, 3D if the
  spline is 3D.

* **use_fista** - Use FISTA deconvolution for peak finding. If this is not set then the
  analysis will be done using a matched filter for peak finding. This is much faster, but
  possibly less accurate at higher densities.

Spliner (EMCCD)
~~~~~~~~~~~~~~~

* **camera_gain** - Conversion factor to go from camera ADU to photo-electrons. Units are ADU/e-,
  so the camera ADU values will be divided by this number to convert to photo-electrons (e-).

* **camera_offset** - This what the camera reads with the shutter closed.

Spliner (sCMOS)
~~~~~~~~~~~~~~~
* **camera_calibration** - This file contains the sCMOS calibration data for the region of
  the camera that the movie comes from. It consists of 4 numpy arrays, [offset, variance, gain,
  relative QE], each of which is the same size as a frame of the movie that is to be analyzed.
  This can be generated for a camera using camera_calibration.py and (if it needs
  to be resliced), reslice_calibration.py.
        
Spliner Standard
~~~~~~~~~~~~~~~~

* **z_value** - Z value(s) in microns at which we will perform convolution with the PSF for
  the purposes of peak finding. If this is not specified the default value is
  z = [0.0]. These are also the starting z values for fitting.

  If your PSF is not degenerate* in Z then it could be helpful to try multiple z
  starting values. However most common 3D PSFs such as astigmatism do not meet
  this criteria. The most commonly used PSF that does meet this criteria is the
  double-helix PSF.

  .. note:: By degenerate I mean that the PSF at one z value can be modeled (with reasonable
	    accuracy) by summing several PSFs with a different z value. For example, most
	    astigmatic PSFs z != 0 can be modeled by summing several z = 0 PSFs with
	    variable x,y positions.

Spliner DECON
~~~~~~~~~~~~~

Spliner using compressed sensing deconvolution for peak finding.

  .. note:: You would typically only do 1 iteration of peak finding and fitting in this case.

* **background_estimator** - Method to use for background estimation, either 'RollingBall' or
  'Wavelet'.

* **decon_method** - Use a compressed sensing deconvolution method for peak finding. If this is
  not specified then peaks are identified by convolving the image with the PSF at one or more
  z values (Spliner Standard). Possible values are 'FISTA' and 'ADMM'

Parameters for 'ADMM' CS deconvolution.

* **admm_iterations** - Iterations of ADMM deconvolution to perform. The larger this value
  is the sharper the peaks will be.

* **admm_lambda** - ADMM lambda value. Larger values will increase the sparsity of the
  deconvolved image.
  
* **admm_number_z** - The number of z-planes to use in the deconvolution. More planes will
  give higher accuracy at the expense of running time, but see the note about z_value in
  spliner standard section as that also applies here.

* **admm_rho** - ADMM rho parameter. A value like 0.1 seems to work well.
  
* **admm_threshold** - Local maxima in the ADMM deconvolved image with values larger than
  this will input into the fitter as localizations to be fit. This number should be roughly
  the minimum peak height that would be considered real times the integral of a peak of this height.

Parameters for 'FISTA' CS deconvolution.

* **fista_iterations** - Iterations of FISTA deconvolution to perform. The larger this value
  is the sharper the peaks will be.

* **fista_lambda** - FISTA lambda value. Larger values will increase the sparsity of the
  deconvolved image.
  
* **fista_number_z** - The number of z-planes to use in the deconvolution. More planes will
  give higher accuracy at the expense of running time, but see the note about z_value in
  spliner standard section as that also applies here.

* **fista_threshold** - Local maxima in the FISTA deconvolved image with values larger than
  this will input into the fitter as localizations to be fit. This number should be roughly
  the minimum peak height that would be considered real times the integral of a peak of this height.

* **fista_timestep** - FISTA timestep. Larger values will cause FISTA to converge faster,
  but if the value is too large FISTA will rapidly diverge.
  
Parameters for 'RollingBall' background removal.

* **rb_radius** - Radius of the rolling ball in pixels.

* **rb_sigma** - Sigma in pixels of the gaussian smoothing to apply to the background
  estimate after the rolling ball step.

Parameters for 'Wavelet' background removal.
            
* **wbgr_iterations** - The number of iterations of background estimation and foreground
  replacement to perform (see the Galloway paper), usually something like 2.

* **wbgr_threshold** - This is the difference between the current estimate and the signal
  at which the signal we be considered "foreground". This should probably be something like 1x
  to 2x the estimated noise in the background.

* **wbgr_wavelet_level** - How many levels of wavelet decomposition to perform. The
  larger the number the less response to local changes in the background, usually something
  like 2.

Multiplane
-----------

This fitter works with any of 3 PSF models (1) the measured PSFs, (2) the pupil functions
or (3) the 3D cubic splines. However you can not mix and match, the models for each
channel/plane must all be of the same type.

* **channelX_cal** - (X = 0-7) The sCMOS camera calibration file for plane X.

* **channelX_ext** - (X = 0-7) The movie file extension for the movie for plane X. The
  analysis works best with a naming scheme like movie_01_c1.tif, movie_01_c2.tif, ...

* **channelX_offset** - (X = 0-7) This parameter allows you to compensate for the
  problem that their might be frame number offsets between the movies from different
  cameras due to synchronization issues.

* **mapping** - The file that contains the transforms for mapping points from one plane
  to another plane.

* **psfX** - (X = 0-7) The PSF files to use for fitting.

* **pupildnX** - (X = 0-7) The pupil function files to use for fitting.

* **splineX** - (X = 0-7) The spline files to use for fitting. These are always 3D splines.
	    
* **weights** - This file contains information about how to weight the per channel/plane
  localization parameters (i.e. x, y, z, etc..) to get the most accurate average value.
  
* **z_value** - Initial z values to consider as starting points for localization z locations.
  Values are in microns.

Pupil Function
--------------

* **pupil_function** - This is the pupil function file to use for fitting.
  
* **z_value** - Initial z values to consider as starting points for localization z locations.
  Values are in microns.
  
Pupil Function (EMCCD)
~~~~~~~~~~~~~~~~~~~~~~

* **camera_gain** - Conversion factor to go from camera ADU to photo-electrons. Units are ADU/e-,
  so the camera ADU values will be divided by this number to convert to photo-electrons (e-).

* **camera_offset** - This what the camera reads with the shutter closed.

Pupil Function (sCMOS)
~~~~~~~~~~~~~~~~~~~~~~
* **camera_calibration** - This file contains the sCMOS calibration data for the region of
  the camera that the movie comes from. It consists of 4 numpy arrays, [offset, variance, gain,
  relative QE], each of which is the same size as a frame of the movie that is to be analyzed.
  This can be generated for a camera using camera_calibration.py and (if it needs
  to be resliced), reslice_calibration.py.
  
PSF FFT
-------

* **psf** - This is the psf file to use for fitting.
  
* **z_value** - Initial z values to consider as starting points for localization z locations.
  Values are in microns.
  
PSF FFT (EMCCD)
~~~~~~~~~~~~~~~

* **camera_gain** - Conversion factor to go from camera ADU to photo-electrons. Units are ADU/e-,
  so the camera ADU values will be divided by this number to convert to photo-electrons (e-).

* **camera_offset** - This what the camera reads with the shutter closed.

PSF FFT (sCMOS)
~~~~~~~~~~~~~~~
* **camera_calibration** - This file contains the sCMOS calibration data for the region of
  the camera that the movie comes from. It consists of 4 numpy arrays, [offset, variance, gain,
  relative QE], each of which is the same size as a frame of the movie that is to be analyzed.
  This can be generated for a camera using camera_calibration.py and (if it needs
  to be resliced), reslice_calibration.py.
	  
L1H
---

* **a_matrix** - A file containing the A matrix to use.

* **epsilon** - Epsilon, `Zhu et al <http://dx.doi.org/doi:10.1038/nmeth.1978>`_ suggest 1.5 for
  poisson simulated data, 2.1 for EMCCD data.
