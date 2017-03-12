Analysis Parameters
===================

These **.xml** files describe how to do the analysis. The parsing
of these files is handled by ``storm_analysis.sa_library.parameters.py``.
Example setting files can be found in ``storm_analysis/test/data``.

General
-------

These parameters are common to all of the analysis programs.

* **append_metadata** - Save the analysis parameters as XML at the end of the
  localization binary file (the list.bin file). If this not specified then
  the analysis parameters will not be appended. ``0`` = No.

*Analysis parameters.*

* **baseline** - This is what the camera reads with the shutter closed.

* **max_frame** - The frame to stop analysis on, ``-1`` = analyze to the end of the film.

* **max_z** - Maximum z value for z fitting, specified in um.
  
* **min_z** - Minimum z value for z fitting, specified in um.

* **orientation** - CCD orientation, generally you should use "normal", but if
  you want to compare the analysis with older versions of Insight3 you'll sometimes
  find that "inverted" works best.

* **peak_locations** - This is for is you already know where your want fitting to
  happen, as for example in a bead calibration movie and you just want to use the
  approximate locations as inputs for fitting.

  peak_locations is a text file with the peak x, y, height and background
  values as white spaced columns (x and y positions are in pixels as
  determined using **visualizer**). ::
  
    1.0 2.0 1000.0 100.0
    10.0 5.0 2000.0 200.0
    ...
  
* **pixel_size** - CCD pixel size (in nm).

* **sigma** - This is the estimated sigma of the PSF in pixels, it serves several
  purposes.

  (1) It is used in most of the analysis approaches as a measure of the
      peak to peak distance at which peak fits do not substantially
      effect each other.

  (2) In most of the analysis, if two peaks are closer than this distance
      then the dimmer one will be discarded.
  
  (3) In **3D-DAOSTORM** and **sCMOS** analysis it is also used as the initial guess
      for the peak sigma.

* **start_frame** - The frame to start analysis on, ``-1`` = start at the beginning of the film.

* **static_background_estimate** - If this is set, and set to a number greater than 0,
  then the analysis will estimate the background by using the average over this number of
  frames.

  If this is not set, or set to 0, the background is estimated separately for each frame.

* **x_start** - X start of the analysis AOI, leave unset to start at the edge of the image.

* **x_stop** - X end of the analysis AOI, leave unset to end at the edge of the image.

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

* **z_correction** - Do z drift correction, ``0`` = No.

3D-DAOSTORM and sCMOS
---------------------

These parameters are common to 3D-DAOSTORM and sCMOS analysis.

* **cutoff** - Z fit cutoff (used when z is calculated later from the fit width in x and y.

* **do_zfit** - Do z fitting (or not), only relevant for "3d" fitting (see "model" parameter).

* **find_max_radius** - To be a peak it must be the maximum value within this radius (in pixels).

* **iterations** - Maximum number of iterations for new peak finding. Usually you'd use
  something like 20 to try and separate neighboring peaks. However if you are analyzing
  beads or some other bright and sparse target 1 may be a better choice as it will suppress
  the spurious splitting of a single peak into multiple peaks.

* **model** - Model is one of "2dfixed", "2d", "3d", or "Z". ::

    2dfixed - fixed sigma 2d gaussian fitting.
    2d - variable sigma 2d gaussian fitting.
    3d - x, y sigma are independently variable, z
         will be fit after peak fitting.
    Z - x, y sigma depend on z, z is fit as part
         of peak fitting.
              
* **wx vs z parameters** - These are used for determining the localization Z position
  based on its in width in x and y (astigmatism imaging). See Huang, Science 2008 for
  a more detailed explanation.
            
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

* **z_step** - The z step size for finding the optimal z value when using the 3d model. If
  this is not specified it defaults to 1.0nm.

3D-DAOSTORM
-----------

* **filter_sigma** - Gaussian filter sigma, this is the sigma of a 2D gaussian to convolve the
  data with prior to peak indentification. When your data has a low SNR this can help for peak
  finding. For optimal sensitivity it should be the same as the expected sigma for your peaks.

  Also, you will need to adjust your threshold parameter as the threshold is now used for
  peak finding in the convolved image and not the original image.
      
  If you set it to zero (or leave it out) then this will not be performed, which can
  make the analysis faster.
	  
* **sigma** - Initial guess for sigma, this is in units of pixels. If you are using the
  2dfixed model then it needs to be pretty close to the correct value. For 2d it should be
  close, probably within 50% or so of the average peak sigma or the fitting might fail
  to converge on many peaks. 3d is similar to 2d. It should not effect fitting for Z
  the model.

  Also see the description of this parameter in *General analysis parameters* section.

* **threshold** - Threshold for a maximum to considered a peak. This is the
  single most important parameter for setting what does and what does not
  get detected in a movie.

  Usually this is the same as the minimum height parameter for peak finding in
  Insight3 (but see note above if you are also using "filter sigma"). You should
  use a number roughly equal to the value of the brightest pixel (minus the CCD
  baseline) in the dimmest peak that you want to detect. If this is too low more
  background will be detected. If it is too high more peaks will be missed.

sCMOS
-----

* **camera_calibration** - This file contains the sCMOS calibration data for the region of
  the camera that the movie comes from. It consists of 3 numpy arrays, [offset, variance, gain],
  each of which is the same size as a frame of the movie that is to be analyzed.
  This can be generated for a camera using camera_calibration.py and (if it needs
  to be resliced), reslice_calibration.py.

* **sigma** - Initial guess for sigma, this is in units of pixels.
  
  For sCMOS analysis it is used as the sigma psf in image segmentation
  (see Section 3.1 of the supplementary material of:
  
  "Video-rate nanoscopy using sCMOS camera-specific single-molecule localization algorithms"
  Huang et al, Nature Methods, 2013.

  It also used to initialize fitting. If you are using the 2dfixed model then it
  needs to be pretty close to the correct value. For 2d it should be close, probably 
  within 50% or so of the average peak sigma or the fitting might fail to converge
  on many peaks. 3d is similar to 2d. It should not effect fitting for Z the model.

* **threshold** - Threshold for a maximum to considered a peak. This has the same meaning
  as for 3D-DAOSTORM, except that it is applied to the convolved image, so you will likely
  need to use a smaller value.
          

Spliner
-------

* **spline** - This is the spline file to use for fitting. Based on the spline the analysis
  will decide whether to do 2D or 3D spline fitting, 2D if the spline is 2D, 3D if the
  spline is 3D.

* **use_fista** - Use FISTA deconvolution for peak finding. If this is not set then the
  analysis will be done using a matched filter for peak finding. This is much faster, but
  possibly less accurate at higher densities.

Spliner standard
~~~~~~~~~~~~~~~~

* **find_max_radius** - To be a peak it must be the maximum value within this radius (in pixels).

* **iterations** - Maximum number of iterations for new peak finding. This has the same
  meaning as for 3D-DAOSTORM and sCMOS analysis.

* **threshold** - Threshold for a maximum to considered a peak. This has the same meaning
  as for 3D-DAOSTORM, except that it is applied to the image convolved with the PSF, so
  you will likely need to use a smaller value.

* **z_value** - Z value(s) in nanometers at which we will perform convolution with the PSF for
  the purposes of peak finding. If this is not specified the default value is
  z = [0.0]. These are also the starting z values for fitting.

  If your PSF is not degenerate* in Z then it could be helpful to try multiple z
  starting values. However most common 3D PSFs such as astigmatism do not meet
  this criteria. The only one that does meet this criteria that is in (sort of)
  common use is the double-helix PSF.

  .. note:: By degenerate I mean that the PSF at one z value can be modeled (with reasonable
	    accuracy) by summing several PSFs with a different z value. For example, most
	    astigmatic PSFs z != 0 can be modeled by summing several z = 0 PSFs with
	    variable x,y positions.

Spliner FISTA
~~~~~~~~~~~~~

FISTA peak finding.

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

* **fista_upsample** - The amount of upsampling to use before FISTA deconvolution. Larger values
  should allow the separation of closer peaks at the expense of running time and (probably)
  speed convergence.

Peak fitting.

* **threshold** - This is basically the same as the minimum height parameter for peak
  finding in Insight3. You should use a number roughly equal to the value of the
  brightest pixel (minus the CCD baseline) in the dimmest peak that you want to
  keep. To some extent this is redundant with the FISTA threshold parameter.
  If you set it much lower than the equivalent FISTA value then it won't make
  much difference. If you set it higher it can remove some of the noise peaks
  that make it through the FISTA step.

* **sigma** - If there are two peaks closer than this value after fitting the dimmer
  one will be removed. Units are in pixels.
          
Rolling Ball background removal. If these are set then this mode of background
estimation will be used (instead of the wavelet based approach below).

* **rb_radius** - Radius of the rolling ball in pixels.

* **rb_sigma** - Sigma in pixels of the gaussian smoothing to apply to the background
  estimate after the rolling ball step.

Wavelet background removal.
            
* **wbgr_iterations** - The number of iterations of background estimation and foreground
  replacement to perform (see the Galloway paper), usually something like 2.

* **wbgr_threshold** - This is the difference between the current estimate and the signal
  at which the signal we be considered "foreground". This should probably be something like 1x
  to 2x the estimated noise in the background.

* **wbgr_wavelet_level** - How many levels of wavelet decomposition to perform. The
  larger the number the less response to local changes in the background, usually something
  like 2.

L1H
---

* **a_matrix** - A file containing the A matrix to use.

* **epsilon** - Epsilon, in Bo's paper he suggested 1.5 for poisson simulated data,
  2.1 for EMCCD data.

