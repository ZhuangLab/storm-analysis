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

* **cutoff** - Z fit cutoff (used when z is calculated later from the fit width in x and y.

* **do_zfit** - Do z fitting (or not), only relevant for "3d" fitting (see "model" parameter).

* **filter_sigma** - Gaussian filter sigma, this is the sigma of a 2D gaussian to convolve the
  data with prior to peak indentification. When your data has a low SNR this can help for peak
  finding. For optimal sensitivity it should be the same as the expected sigma for your peaks.

  Also, you will need to adjust your threshold parameter as the threshold is now used for
  peak finding in the convolved image and not the original image.
      
  If you set it to zero (or leave it out) then this will not be performed, which can
  make the analysis faster.

  .. note:: This is not relevant for sCMOS analysis.

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
  

            
            # wx vs z parameters
            #
            # See Huang, Science 2008 for a more detailed explanation.
            #
            "wx_wo" : ["float", None],
            "wx_c" : ["float", None],
            "wx_d" : ["float", None],
            "wxA" : ["float", None],
            "wxB" : ["float", None],
            "wxC" : ["float", None],
            "wxD" : ["float", None],

            # wy vs z parameters.
            "wy_wo" : ["float", None],
            "wy_c" : ["float", None],
            "wy_d" : ["float", None],
            "wyA" : ["float", None],
            "wyB" : ["float", None],
            "wyC" : ["float", None],
            "wyD" : ["float", None],

            # The starting z value for fitting. If this is not specified it defaults to 0.0.
            "z_value" : ["float", None],

            # The z step size for finding the optimal z value when using the 3d model. If
            # this is not specified it defaults to 1.0.
            "z_step" : ["float", None],

L1H
---

sCMOS
-----

Spliner
-------

Spliner standard
~~~~~~~~~~~~~~~~

Spliner FISTA
~~~~~~~~~~~~~
