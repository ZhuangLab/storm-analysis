#!/usr/bin/env python
"""
This contains some of what is common to all of the peak finding algorithms.
 

Note on handling RQE differences: 

1) The images are corrected at loading for RQE differences, which is done 
by dividing the image be the measured camera RQE.

2) The peak finding classes use this corrected image. I think this is not 
quite correct from a statistical perspective, but is simple and works well 
for the small RQE differences (~5%) that are usually encountered. "Works
well" means that the finder is noticeably biased against localizations in
areas with lower than average RQE.

3) The peak fitting (C library) undoes the RQE correction, it multiplies
the image it receives by the RQE correction. This way the fitting is correct,
at least in a statistical sense.

Hazen 06/19
"""
import numpy
import os
import tifffile

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.draw_gaussians_c as dg

#
# Functions.
#
def gaussianPSF(shape, sigma):
    """
    Return a normalized 2D Gaussian, usually used for creating MatchedFilter objects.
    """
    psf = dg.drawGaussiansXY(shape,
                             numpy.array([0.5*shape[0]]),
                             numpy.array([0.5*shape[1]]),
                             sigma = sigma)
    return psf/numpy.sum(psf)

def getPeakLocations(peak_filename, margin, pixel_size, sigma):
    """
    This is for if you already know where your want fitting to happen, as
    for example in a bead calibration movie and you just want to use the
    approximate locations as inputs for fitting.

    There are two choices for peak_locations file format:

    1. A text file with the peak x, y, height and background values as 
       white spaced columns (x and y positions are in pixels as determined 
       using visualizer).

       1.0 2.0 1000.0 100.0
       10.0 5.0 2000.0 200.0
       ...

    2. An HDF5 format localization file. This is treated in a similar
       fashion to the text file in that all of the locations are loaded.
       If the fields 'xsigma' or 'ysigma' exist they will be used for
       the initial X/Y sigma values of the localization.
    """
    if os.path.exists(peak_filename):
        print("Using peak starting locations specified in", peak_filename)
    elif os.path.exists(os.path.basename(peak_filename)):
        peak_filename = os.path.basename(peak_filename)
        print("Using peak starting locations specified in", peak_filename)

    # Check if the file is a storm-analysis HDF5 file.
    #
    if saH5Py.isSAHDF5(peak_filename):
        peak_locations_type = "hdf5"
        peak_locations = saH5Py.loadLocalizations(peak_filename)
        if not "ysigma" in peak_locations:
            if not "xsigma" in peak_locations:
                peak_locations["xsigma"] = numpy.ones(peak_locations["x"].size) * sigma
            peak_locations["ysigma"] = peak_locations["xsigma"].copy()
            
    else:
        peak_locations_type = "text"

        # Load peak x,y locations.
        peak_locs = numpy.loadtxt(peak_filename, ndmin = 2)

        # Create peak dictionary.
        peak_locations = {"background" : peak_locs[:,3],
                          "height" : peak_locs[:,2],
                          "x" : peak_locs[:,0],
                          "y" : peak_locs[:,1]}

        peak_locations["xsigma"] = numpy.ones(peak_locations["x"].size) * sigma
        peak_locations["ysigma"] = numpy.ones(peak_locations["x"].size) * sigma
        peak_locations["z"] = numpy.zeros(peak_locations["x"].size)

    # Adjust positions for finding/fitting margin.
    peak_locations["x"] += margin
    peak_locations["y"] += margin

    print("Loaded", peak_locations["x"].size, "peak locations")
    #
    # We return is_text as the caller might want to do different things if
    # the file is text, like initialize the Z value.
    #
    return [peak_locations, peak_locations_type]

def padArray(ori_array, pad_size):
    """
    Pads out an array to a large size.

    ori_array - A 2D numpy array.
    pad_size - The number of elements to add to each of the "sides" of the array.
    
    The padded 2D numpy array.
    """
    if (pad_size > 0):
        [x_size, y_size] = ori_array.shape
        lg_array = numpy.ones((x_size+2*pad_size,y_size+2*pad_size))
        lg_array[pad_size:(x_size+pad_size),pad_size:(y_size+pad_size)] = ori_array.astype(numpy.float64)
        lg_array[0:pad_size,:] = numpy.flipud(lg_array[pad_size:2*pad_size,:])
        lg_array[(x_size+pad_size):(x_size+2*pad_size),:] = numpy.flipud(lg_array[x_size:(x_size+pad_size),:])
        lg_array[:,0:pad_size] = numpy.fliplr(lg_array[:,pad_size:2*pad_size])
        lg_array[:,(y_size+pad_size):(y_size+2*pad_size)] = numpy.fliplr(lg_array[:,y_size:(y_size+pad_size)])
        return lg_array
    
    else:
        return ori_array

def peakMask(shape, parameters, margin):
    """
    Return the array that is used to mask the image to reduce the ROI
    where peaks can be found.
    """
    peak_mask = numpy.ones(shape)

    # Check for circular AOI.
    if parameters.hasAttr("x_center"):
        assert parameters.hasAttr("y_center"), "Y center must be specified."
        assert parameters.hasAttr("aoi_radius"), "AOI radius must be specified."

        rr = parameters.getAttr("aoi_radius")
        xc = parameters.getAttr("x_center") + margin
        yc = parameters.getAttr("y_center") + margin

        rr = rr*rr
        xr = numpy.arange(peak_mask.shape[0]) - xc
        yr = numpy.arange(peak_mask.shape[0]) - yc
        xv, yv = numpy.meshgrid(xr, yr)
        peak_mask[((xv*xv + yv*yv) > rr)] = 0

    # This catches 'y_center' without 'x_center'.
    elif parameters.hasAttr("y_center"):
        assert parameters.hasAttr("x_center"), "X center must be specified."

    # This catches 'aoi_radius' without 'x_center'.
    elif parameters.hasAttr("aoi_radius"):
        assert parameters.hasAttr("x_center"), "X center must be specified."
        
    # Check for square AOI
    else:
        if parameters.hasAttr("x_start"):
            peak_mask[:,0:parameters.getAttr("x_start")+margin] = 0.0
        if parameters.hasAttr("x_stop"):
            peak_mask[:,parameters.getAttr("x_stop")+margin:-1] = 0.0
        if parameters.hasAttr("y_start"):
            peak_mask[0:parameters.getAttr("y_start")+margin,:] = 0.0
        if parameters.hasAttr("y_stop"):
            peak_mask[parameters.getAttr("y_stop")+margin:-1,:] = 0.0
        
    return peak_mask


#
# Classes.
#
def FittingException(Exception):
    pass


class PeakFinder(object):
    """
    Base class for peak finding. This handles identification of peaks in an image.

    If you want to modify this with a custom peak finder or an alternative
    way to estimate the background, the recommended approach is to sub-class this
    class and then modify backgroundEstimator(), newImage() and peakFinder().
    """
    def __init__(self, parameters = None, **kwds):
        """
        This is called once at the start of analysis to initialize the
        parameters that will be used for peak fitting.
 
        parameters - A parameters object.
        """
        super(PeakFinder, self).__init__(**kwds)
        
        # Initialized from parameters.
        self.find_max_radius = parameters.getAttr("find_max_radius")     # Radius (in pixels) over which the maxima is maximal.
        self.iterations = parameters.getAttr("iterations")               # Maximum number of cycles of peak finding, fitting and subtraction to perform.
        self.sigma = parameters.getAttr("sigma")                         # Peak sigma (in pixels).
        self.threshold = parameters.getAttr("threshold")                 # Peak minimum threshold in units of sigma (as in "3 sigma effect").
        
        # Other member variables.
        self.background = None                                           # Current estimate of the image background.
        self.bg_filter = None                                            # Background MatchedFilter object.
        self.camera_variance = None                                      # Camera variance, only relevant for a sCMOS camera.
        self.check_mode = False                                          # Run in diagnostic mode. Only useful for debugging.
        self.image = None                                                # The original image.
        self.margin = None                                               # Size of the unanalyzed "edge" around the image.
        self.mfinder = None                                              # The maxima finder.
        self.parameters = parameters                                     # Keep access to the parameters object.
        self.peak_locations = None                                       # Initial peak locations, as explained below.
        self.peak_locations_type = None                                  # Initial peak locations type.
        self.peak_mask = None                                            # Mask for limiting peak identification to a particular AOI.
        self.pf_iterations = 0                                           # Keep track of the total number of iterations that were performed.
        
        # Print warning about check mode
        if self.check_mode:
            print("Warning! Running in check mode!")

        # Only do one cycle of peak finding as we'll always return the same locations.
        if parameters.hasAttr("peak_locations"):
            if (self.iterations != 1):
                print("WARNING: setting number of iterations to 1!")
                self.iterations = 1            
        
    def backgroundEstimator(self, image):
        """
        This method does the actual background estimation. It is just a simple
        low pass filter.

        Override this if you want to change how the background is estimated.

        FIXME: Convolution should be weighted by the camera variance if it is 
               not uniform?
        """
        return self.bg_filter.convolve(image)

    def cleanUp(self):
        print("  ", self.pf_iterations, "peak finding iterations.")
        print()
        if self.bg_filter is not None:
            self.bg_filter.cleanup()

    def estimateBackground(self, fit_peaks_image, bg_estimate):
        """
        Estimate the background for the image.
        
        fit_peaks_image - The current best fit image.
        bg_estimate - An estimate of the background.

        Returns the current background estimate.
        """        
        # If we are provided with an estimate of the background then just use it.
        if bg_estimate is not None:
            self.background = bg_estimate

        # Otherwise make our own estimate.
        else:
            image = self.image - fit_peaks_image
            self.background = self.backgroundEstimator(image)

        if self.check_mode:
            with tifffile.TiffWriter("bg_estimate.tif") as tf:
                tf.save(self.background.astype(numpy.float32))

        return self.background

    def findPeaks(self, fit_peaks_image):
        """
        Finds the peaks in the image.
   
        fit_peaks_image - The current fit image.
    
        Return [new peaks, new peak type, done]
        """
        self.pf_iterations += 1
        
        # Use pre-specified peak locations if available, e.g. bead calibration.
        if self.peak_locations is not None:
            return [self.peak_locations, self.peak_locations_type, True]
            
        # Otherwise, identify local maxima in the image.
        new_peaks = self.peakFinder(fit_peaks_image)

        # Update new peak identification threshold (if necessary).
        # Also, while threshold is greater than min_threshold we
        # are automatically not done.
        if (self.cur_threshold > self.threshold):
            self.cur_threshold -= 1.0
            return [new_peaks, "finder", False]

        # If we did not find any new peaks then we may be done.
        if (new_peaks["x"].size == 0):
            return [new_peaks, "finder", True]
        else:
            return [new_peaks, "finder", False]

    def newImage(self, new_image):
        """
        This is called once at the start of the analysis of a new image.
        
        new_image - A 2D numpy array.
        """
        # Make a copy of the starting image.
        #
        # FIXME: Is this copy necessary? We're not doing this in multiplane.
        #
        self.image = numpy.copy(new_image)

        # Initialize new peak minimum threshold. If we are doing more
        # than one iteration we start a bit higher and come down to
        # the specified threshold.
        if(self.iterations>4):
            self.cur_threshold = self.threshold + 4.0
        else:
            self.cur_threshold = self.threshold + float(self.iterations)

        # Create mask to limit peak finding to a user defined sub-region of the image.
        if self.peak_mask is None:
            self.peak_mask = peakMask(new_image.shape, self.parameters, self.margin)

        # Create filter objects if necessary.
        if self.bg_filter is None:

            # Create matched filter for background.
            bg_psf = gaussianPSF(new_image.shape, self.parameters.getAttr("background_sigma"))
            self.bg_filter = matchedFilterC.MatchedFilter(bg_psf,
                                                          fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                          memoize = True,
                                                          max_diff = 1.0e-3)

        # Reset maxima finder.
        self.mfinder.resetTaken()

    def padArray(self, np_array):
        """
        Return a version of array padded to the correct size.
        """
        return padArray(np_array, self.margin)
    
    def peakFinder(self, fit_peaks_image):
        """
        Sub-classes must provide this method.
        """
        raise FittingException("Finder had no peakFinder() method.")

    def setVariance(self, camera_variance):
        """
        Set the camera variance, usually used in sCMOS analysis.
        """
        self.camera_variance = self.padArray(camera_variance)
        return self.camera_variance
            
    
class PeakFinderGaussian(PeakFinder):
    """
    This is the peak finder for 3D-DAOSTORM and sCMOS, it handles Gaussian shaped peaks.
    """
    def __init__(self, parameters = None, **kwds):
        """
        This is called once at the start of analysis to initialize the
        parameters that will be used for peak fitting.
 
        parameters - A parameters object.
        """
        kwds["parameters"] = parameters
        super(PeakFinderGaussian, self).__init__(**kwds)

        # Figure out what margin and ROI to use.
        if (self.parameters.getAttr("roi_size", -1) != -1):
            self.roi_size = parameters.getAttr("roi_size")
        else:

            # Calculate roi size based on sigma.
            self.roi_size = int(8.0 * self.sigma)

            # Make it even larger for variable width fitters.
            if(parameters.getAttr("model") != "2dfixed"):
                self.roi_size = int(1.5 * self.roi_size)
        self.margin = int(self.roi_size/2 + 2)
        
        # Initialized from parameters.
        self.z_value = self.parameters.getAttr("z_value", 0.0)           # The starting z value to use for peak fitting.
        
        # Other member variables.
        self.fg_mfilter = None                                           # Foreground MatchedFilter object (may be None).
        self.fg_vfilter = None                                           # Foreground variance MatchedFilter object, will
                                                                         # be none if self.fg_mfilter is None.

        # Configure maxima finder.
        #
        self.mfinder = iaUtilsC.MaximaFinder(margin = self.margin,
                                             radius = self.find_max_radius,
                                             threshold = self.threshold,
                                             z_values = [self.z_value])
        
        # Load peak locations if specified.
        #
        # FIXME: The starting z value is always 0.0. Not sure why we don't use
        #        self.z_value for this. Though I guess it would only really be
        #        relevant for the 'Z' fitting model.
        #
        if parameters.hasAttr("peak_locations"):
            [self.peak_locations, self.peak_locations_type] = getPeakLocations(parameters.getAttr("peak_locations"),
                                                                               self.margin,
                                                                               parameters.getAttr("pixel_size"),
                                                                               self.sigma)

    def cleanUp(self):
        super(PeakFinderGaussian, self).cleanUp()
        if self.fg_mfilter is not None:
            self.fg_mfilter.cleanup()
            self.fg_vfilter.cleanup()

    def getROISize(self):
        return self.roi_size

    def newImage(self, new_image):
        super(PeakFinderGaussian, self).newImage(new_image)

        #
        # Create matched filter for foreground as well as a matched filter
        # for calculating the expected variance of the background if it was
        # smoothed on the same scale as the foreground.
        #
        if (self.fg_mfilter is None) and self.parameters.hasAttr("foreground_sigma"):
            if (self.parameters.getAttr("foreground_sigma") > 0.0):
                fg_psf = gaussianPSF(new_image.shape, self.parameters.getAttr("foreground_sigma"))
                self.fg_mfilter = matchedFilterC.MatchedFilter(fg_psf,
                                                               fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                               memoize = True,
                                                               max_diff = 1.0e-3)
                self.fg_vfilter = matchedFilterC.MatchedFilter(fg_psf * fg_psf,
                                                               fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                               memoize = True,
                                                               max_diff = 1.0e-3)

    def peakFinder(self, fit_peaks_image):
        """
        This method does the actual peak finding.
        """
        # Calculate background variance.
        #
        # Note the assumption here that we are working in units of photo-electrons
        # so Poisson statistics applies, variance = mean.
        #
        bg_var = self.background + fit_peaks_image
        
        # Add camera variance if set.
        if self.camera_variance is not None:
            bg_var += self.camera_variance

        # Calculate weighted variance if the image is being smoothed.
        if self.fg_vfilter is not None:
            bg_var = self.fg_vfilter.convolve(bg_var)

        if self.check_mode:
            with tifffile.TiffWriter("variances.tif") as tf:
                tf.save(bg_var.astype(numpy.float32))
            
        # Remove problematic values.
        #
        mask = (bg_var <= 0.1)
        if (numpy.sum(mask) > 0):
            if self.check_mode:
                print("Warning! small and/or negative values detected in background variance!")
            bg_var[mask] = 0.1
                
        # Convert to standard deviation.
        bg_std = numpy.sqrt(bg_var)

        # Calculate foreground.
        foreground = self.image - self.background - fit_peaks_image

        # Calculate smoothed image if we have a foreground filter.
        if self.fg_mfilter is not None:
            foreground = self.fg_mfilter.convolve(foreground)

        if self.check_mode:
            with tifffile.TiffWriter("foreground.tif") as tf:
                tf.save(foreground.astype(numpy.float32))
            
        # Calculate foreground in units of signal to noise.
        foreground = foreground/bg_std
                    
        if self.check_mode:
            with tifffile.TiffWriter("fg_bg_ratio.tif") as tf:
                tf.save(foreground.astype(numpy.float32))
        
        # Mask the image so that peaks are only found in the AOI.
        masked_image = foreground * self.peak_mask

        # Identify local maxima in the masked image.
        [x, y, z] = self.mfinder.findMaxima([masked_image])
        return {"x" : x, "y" : y, "z" : z, "sigma" : numpy.ones(x.size)*self.sigma}


class PeakFinderArbitraryPSF(PeakFinder):
    """
    This is the base class for Spliner, Pupilfn and PSFFFT, it handles arbitrary
    PSF shapes possibly with multiple z values.
    """
    def __init__(self, parameters = None, psf_object = None, **kwds):
        kwds["parameters"] = parameters
        super(PeakFinderArbitraryPSF, self).__init__(**kwds)

        self.fg_mfilter = []
        self.fg_mfilter_zval = []
        self.fg_vfilter = []
        self.psf_object = psf_object
        self.z_values = []

        self.margin = psf_object.getMargin()
        
        # Note: self.z_values is the Z position in 'internal' units, i.e. the units
        #       that the PSF generation library uses. For splines for example this
        #       is the spline size in Z.
        #
        # 'z_value' is in units of microns.
        #  self.fg_mfilter_zval is the Z position in nanometers.
        #
        self.fg_mfilter_zval = list(map(lambda x: x * 1.0e+3, parameters.getAttr("z_value", [0.0])))
        for zval in self.fg_mfilter_zval:
            assert self.psf_object.isValidZ(zval)
            self.z_values.append(self.psf_object.getScaledZ(zval))

        # Configure maxima finder.
        #
        self.mfinder = iaUtilsC.MaximaFinder(margin = self.margin,
                                             radius = self.find_max_radius,
                                             threshold = self.threshold,
                                             z_values = self.z_values)

        if parameters.hasAttr("peak_locations"):
            [self.peak_locations, self.peak_locations_type] = getPeakLocations(parameters.getAttr("peak_locations"),
                                                                               self.margin,
                                                                               parameters.getAttr("pixel_size"),
                                                                               self.sigma)

            # Set initial z value (for text files).
            if (self.peak_locations_type == "text"):
                self.peak_locations["z"][:] = self.z_values[0]

            # Convert z value to PSF units (for HDF5 localization files).
            else:
                # If the HDF5 file does not have any "z" information we use the first
                # value of self.z_values.
                #
                if not "z" in self.peak_locations:
                    self.peak_locations["z"] = numpy.zeros(self.peak_locations["x"].size)
                    self.peak_locations["z"][:] = self.z_values[0]
                self.peak_locations["z"] = self.psf_object.getScaledZ(self.peak_locations["z"])

    def cleanUp(self):
        super(PeakFinderArbitraryPSF, self).cleanUp()
        for i in range(len(self.fg_mfilter)):
            self.fg_mfilter[i].cleanup()
            self.fg_vfilter[i].cleanup()

    def newImage(self, new_image):
        """
        This is called once at the start of the analysis of a new image.
        
        new_image - A 2D numpy array.
        """
        super(PeakFinderArbitraryPSF, self).newImage(new_image)
    
        #
        # If does not already exist, create filter objects from
        # the PSF at different z values.
        #
        if (len(self.fg_mfilter) == 0):
            for zval in self.fg_mfilter_zval:
                psf = self.psf_object.getPSF(zval,
                                             shape = new_image.shape,
                                             normalize = False)
                psf_norm = psf/numpy.sum(psf)
                fg_mfilter = matchedFilterC.MatchedFilter(psf_norm,
                                                          fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                          memoize = True,
                                                          max_diff = 1.0e-3)
                self.fg_mfilter.append(fg_mfilter)
                self.fg_vfilter.append(matchedFilterC.MatchedFilter(psf_norm * psf_norm,
                                                                    fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                                    memoize = True,
                                                                    max_diff = 1.0e-3))

                # Save a picture of the PSF for debugging purposes.
                if self.check_mode:
                    print("psf max", numpy.max(psf))
                    filename = "psf_{0:.3f}.tif".format(zval)
                    tifffile.imsave(filename, psf.astype(numpy.float32))
                        
    def peakFinder(self, fit_peaks_image):
        """
        This method does the actual peak finding.
        """
        if self.check_mode:
            tifffile.imsave("fit_peaks.tif", fit_peaks_image.astype(numpy.float32))

        # Calculate background variance.
        #
        # Notes:
        #
        # 1. The assumption here is that we are working in units of photo-electrons
        #    so Poisson statistics applies, variance = mean.
        #
        # 2. We use the absolute value of fit_peaks_image as some fitters (such as
        #    Spliner) will sometimes return this array with negative values. This is
        #    still probably not the correct way to handle this.
        #
        bg_var = self.background + numpy.abs(fit_peaks_image)

        # Add camera variance if set.
        if self.camera_variance is not None:
            bg_var += self.camera_variance
            
        #
        # Find peaks in image convolved with the PSF at different z values.
        #
        if self.check_mode:
            bg_tif = tifffile.TiffWriter("background.tif")
            fg_tif = tifffile.TiffWriter("foreground.tif")
            fg_bg_ratio_tif = tifffile.TiffWriter("fg_bg_ratio.tif")

        masked_images = []
        for i in range(len(self.fg_mfilter)):

            # Estimate background variance at this particular z value.
            background = self.fg_vfilter[i].convolve(bg_var)

            # Remove problematic values.
            #
            mask = (background <= 0.1)
            if (numpy.sum(mask) > 0):
                if self.check_mode:
                    print("Warning! small and/or negative values detected in background variance!")
                background[mask] = 0.1
                    
            # Convert to standard deviation.
            bg_std = numpy.sqrt(background)

            # Calculate foreground.
            foreground = self.image - self.background - fit_peaks_image
            foreground = self.fg_mfilter[i].convolve(foreground)

            # Calculate foreground in units of signal to noise.
            fg_bg_ratio = foreground/bg_std
        
            if self.check_mode:
                bg_tif.save(background.astype(numpy.float32))
                fg_tif.save(foreground.astype(numpy.float32))
                fg_bg_ratio_tif.save(fg_bg_ratio.astype(numpy.float32))

            # Mask the image so that peaks are only found in the AOI.
            masked_images.append(fg_bg_ratio * self.peak_mask)

        if self.check_mode:
            bg_tif.close()
            fg_tif.close()
            fg_bg_ratio_tif.close()

        # Identify local maxima in the masked ratio images stack.
        [x, y, z] = self.mfinder.findMaxima(masked_images)
        return {"x" : x, "y" : y, "z" : z}


class PeakFitter(object):
    """
    Base class for peak fitting. This handles refinement of the
    parameters of the peaks that were identified with PeakFinder.

    The actual fitting is done by an the self.mfitter object, this
    is primarily just a wrapper for the self.mfitter object.
    """
    def __init__(self, mfitter = None, parameters = None, **kwds):
        """
        parameters - A (fitting) parameters object.
        """
        super(PeakFitter, self).__init__(**kwds)

        self.no_fitting = parameters.getAttr("no_fitting", 0) # If this is True then we won't do any fitting
                                                              # iterations. This is useful for testing the
                                                              # finder, as well as how accurately we're
                                                              # initializing the peak parameter values.
        self.image = None              # The image for peak fitting.
        self.mfitter = mfitter         # An instance of a sub-class of the MultiFitter class.
        self.minimum_significance = parameters.getAttr("threshold")         # The threshold value is also the minimum peak significance value.
        self.sigma = parameters.getAttr("sigma")                            # Peak sigma (in pixels).
        self.neighborhood = self.sigma * PeakFinderFitter.unconverged_dist  # Radius for marking neighbors as unconverged.

    def cleanUp(self):
        self.mfitter.cleanup()

    def fitPeaks(self, new_peaks, peaks_type):
        """
        Performs a single iteration of peak fitting.
        
        new_peaks - A dictionary containing numpy arrays specifying peak x/y location, etc.
                    to add to the fit.
        peaks_type - The type of the peaks, e.g. 'finder' if they were identified by the 
                     peak finder, or '??' if they were provided by the user.
    
        returns the current fit peaks image.
        """
        # Check if we need to do anything.
        if (new_peaks["x"].size > 0):

            # Update status of current peaks (if any) that are near
            # to the new peaks that are being added.
            #
            if (self.mfitter.getNFit() > 0):
                c_x = self.mfitter.getPeakProperty("x")
                c_y = self.mfitter.getPeakProperty("y")
                status = self.mfitter.getPeakProperty("status")
                new_status = iaUtilsC.runningIfHasNeighbors(status,
                                                            c_x,
                                                            c_y,
                                                            new_peaks["x"],
                                                            new_peaks["y"],
                                                            self.neighborhood)
                self.mfitter.setPeakStatus(new_status)
                
            # Add new peaks.
            self.mfitter.newPeaks(new_peaks, peaks_type)

            # Iterate fitting and remove any error peaks.
            #
            # The assumption is that because error peaks are longer in the
            # fit image we don't have to do additional iterations on the
            # remaining peaks after the error peaks have been removed.
            #
            if not self.no_fitting:
                self.mfitter.doFit()
                self.mfitter.removeErrorPeaks()

            # Remove peaks that are too close to each other and/or that
            # have a low significance score.
            #
            status = self.mfitter.getPeakProperty("status")

            # Identify peaks that are to close based on the somewhat
            # arbitrary criteria of being within 1 sigma.
            #
            # markDimmerPeaks() will update the status array, in particular
            # it will mark the dimmer of two peaks that are too close as ERROR.
            #
            px = self.mfitter.getPeakProperty("x")
            py = self.mfitter.getPeakProperty("y")
            n_proximity = iaUtilsC.markDimmerPeaks(px,
                                                   py,
                                                   self.mfitter.getPeakProperty("height"),
                                                   status,
                                                   self.sigma,
                                                   self.neighborhood)

            # Identify peaks that have a low significance score.
            #
            # markLowSignificancePeaks() will update the status array, in particular
            # it will mark low significance peaks as ERROR.
            #
            n_significance = iaUtilsC.markLowSignificancePeaks(px,
                                                               py,
                                                               self.mfitter.getPeakProperty("significance"),
                                                               status,
                                                               self.minimum_significance,
                                                               self.neighborhood)

            # This does the actual peak removal. We update the peak status in
            # mfitter, then tell mfitter to remove all the ERROR peaks.
            #
            if ((n_proximity + n_significance) > 0):
                self.mfitter.setPeakStatus(status)
                self.mfitter.removeErrorPeaks()
                self.mfitter.incProximityCounter(n_proximity)
                self.mfitter.incSignificanceCounter(n_significance)

            # If we have unconverged peaks, iterate some more.
            if (self.mfitter.getUnconverged() > 0) and (not self.no_fitting):
                self.mfitter.doFit()
                self.mfitter.removeErrorPeaks()

        # Return the current fit image.
        return self.mfitter.getFitImage()

    def getPeakProperty(self, pname):
        return self.mfitter.getPeakProperty(pname)
        
    def newBackground(self, new_background):
        """
        new_background - A new estimate of the image background (2D numpy array).
        """
        self.mfitter.newBackground(new_background)
        
    def newImage(self, new_image):
        """
        new_image - A new image (2D numpy array).

        Note - Addition of the sCMOS calibration term to the image is done
               in the C library.
        """
        if not self.mfitter.isInitialized():
            self.mfitter.initializeC(new_image)
        self.mfitter.newImage(new_image)

    def removeRunningPeaks(self):
        """
        Remove all the unconverged peaks.
        """
        if not self.no_fitting:
            self.mfitter.removeRunningPeaks()
        
    def rescaleZ(self, z):
        """
        Convert from fitting z units to microns.
        """
        return self.mfitter.rescaleZ(z)


class PeakFitterArbitraryPSF(PeakFitter):
    """
    Class for arbitrary PSF based peak fitting.
    """
    def __init__(self, **kwds):
        super(PeakFitterArbitraryPSF, self).__init__(**kwds)

        # Update refitting neighborhood parameter.
        self.neighborhood = int(0.5 * self.mfitter.getSize()) + 1

    
class PeakFinderFitter(object):
    """
    Base class to encapsulate peak finding and fitting. 

    To get an idea of how all the pieces are supposed to go together see:
      3d_daostorm/find_peaks.py
      sCMOS/find_peaks.py
    """
    unconverged_dist = 5.0  # Distance between peaks for marking as unconverged in
                            # pixels, this is multiplied by parameters.sigma. Hopefully
                            # peaks that are > 5 sigma away from each other have little
                            # effect on each other.
    
    def __init__(self, peak_finder = None, peak_fitter = None, properties = None, **kwds):
        """
        peak_finder - A PeakFinder object.
        peak_fitter - A PeakFitter object.
        properties - Which properties to return, e.g. "x", "y", etc.
        """
        super(PeakFinderFitter, self).__init__(**kwds)

        self.peak_finder = peak_finder
        self.peak_fitter = peak_fitter
        self.properties = properties

        # The properties must include at least 'x' and 'y'.
        assert ("x" in self.properties), "'x' is a required property."
        assert ("y" in self.properties), "'y' is a required property."

    def analyzeImage(self, movie_reader):
        """
        Analyze an image and return all the peaks that were found and fit.

        movie_reader - analysis_io.MovieReader object.

        return - found peaks dictionary.
        """
        # Load image (in photo-electrons).
        [image, fit_peaks_image] = self.loadImage(movie_reader)

        # Load background estimate (in photo-electrons).
        bg_estimate = self.loadBackgroundEstimate(movie_reader)

        self.peak_finder.newImage(image)
        self.peak_fitter.newImage(image)

        for i in range(self.peak_finder.iterations):

            # Update background estimate.
            background = self.peak_finder.estimateBackground(fit_peaks_image, bg_estimate)

            # Find new peaks.
            [new_peaks, peaks_type, done] = self.peak_finder.findPeaks(fit_peaks_image)

            # Fit new peaks.
            self.peak_fitter.newBackground(background)
            fit_peaks_image = self.peak_fitter.fitPeaks(new_peaks, peaks_type)
                    
            if done:
                break

        # Remove any peaks that have not converged.
        self.peak_fitter.removeRunningPeaks()

        # Return a dictionary with the requested properties.
        return self.getPeakProperties()

    def cleanUp(self):
        self.peak_finder.cleanUp()
        self.peak_fitter.cleanUp()

    def getPeakProperties(self):
        """
        Create a dictionary with the requested properties.
        """
        peaks = {}
        for pname in self.properties:
            peaks[pname] = self.peak_fitter.getPeakProperty(pname)

            # x,y,z corrections.
            if (pname == "x"):
                peaks[pname] -= float(self.peak_finder.margin)

            elif (pname == "y"):
                peaks[pname] -= float(self.peak_finder.margin)
                
            elif (pname == "z"):
                peaks[pname] = self.peak_fitter.rescaleZ(peaks[pname])

        return peaks
        
    def loadBackgroundEstimate(self, movie_reader):
        bg_estimate = movie_reader.getBackground()
        if bg_estimate is not None:
            bg_estimate = padArray(bg_estimate, self.peak_finder.margin)
            
        return bg_estimate
        
    def loadImage(self, movie_reader):
        image = padArray(movie_reader.getFrame(), self.peak_finder.margin)
        fit_peaks_image = numpy.zeros(image.shape)
        return [image, fit_peaks_image]


class PSFFunction(object):
    """
    This is the base class for handling the PSF for fitters that use
    arbitrary PSFs such as PSFFFT, PupilFN and Spliner. In theory it
    handles all of the details of the PSF, such as how many pixels
    it covers, how much margin to add to the image, how to convert
    Z values, etc..
    """
    def getCPointer(self):
        """
        Returns a pointer to the C library structure that is
        used to describe the PSF.

        pupilData in pupilfn/pupil_function.h for example.
        """
        assert False

    def getMargin(self):
        """
        Return the margin to add to the image for finding/fitting.
        """
        assert False
        
    def getPSF(self, z_value, shape = None, normalize = False):
        """
        Return an image of the PSF at z_value, centered in an 
        array of size shape.
        """
        assert False

    def getScaledZ(self, z_value):
        """
        This expects z_value to be in nanometers.
        """
        assert False
        
    def getSize(self):
        """
        Return the X/Y size in pixels (all the fitters expect the 
        PSF to be square).
        """
        assert False

    def getZMax(self):
        """
        Return maximum z position for the PSF in nanometers.
        """
        return self.zmax        
        
    def getZMin(self):
        """
        Return the minimum z position for the PSF in nanometers.
        """
        return self.zmin

    def isValidZ(self, z):
        """
        Return True if the z value is within the z range covered by the
        PSF (not including the end-points as these are troublesome for spliner).
        """
        if (z <= self.getZMin()):
            return False
        if (z >= self.getZMax()):
            return False
        else:
            return True

    def rescaleZ(self, z_value):
        """
        Convert from fitting units back to *microns* (not nanometers).
        """
        return z_value

#
# The MIT License
#
# Copyright (c) 2014 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
