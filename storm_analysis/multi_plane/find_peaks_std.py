#!/usr/bin/env python
"""
(Standard) multi-plane peak finder and fitter.

Peaks are identified using ideas from this paper:

"SNSMIL, a real-time single molecule identification and localization algorithm 
for super-resolution fluorescence microscopy", Tang et al., Scientific Reports, 2015.

"channel", "frame", "image" and "plane" are used somewhat interchangeably..

FIXME: For multi-plane data with the planes far enough part we should allow
       two peaks to have a similar x,y location if their z locations are
       different enough.

FIXME: Overly complicated? Would it work just as well to transform all the
       input images and sCMOS calibration data first, then do the fitting
       on the transformed image? This would be a lot simpler..

Hazen 05/17
"""
import os
import pickle
import numpy
import tifffile

import storm_analysis.multi_plane.mp_fit_c as mpFitC
import storm_analysis.multi_plane.mp_utilities as mpUtil

import storm_analysis.sa_library.affine_transform_c as affineTransformC
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC

import storm_analysis.psf_fft.psf_fn as psfFn
import storm_analysis.pupilfn.pupil_fn as pupilFn
import storm_analysis.spliner.spline_to_psf as splineToPSF


class MPPeakFinder(fitting.PeakFinder):
    """
    Multi-plane peak finding.

    All locations are relative to plane 0. The number of localizations
    is a multiple of the number of focal planes.

    The expectation is that we are working with an image that has already 
    been corrected for gain and offset, i.e. we are working in units of
    photo-electrons.

    This works with affine transformed versions of the images so that they
    can all be overlaid on top of each other.
    """
    def __init__(self, parameters = None, psf_objects = None, **kwds):
        kwds["parameters"] = parameters
        super(MPPeakFinder, self).__init__(**kwds)

        self.atrans = [None]
        self.backgrounds = []
        self.height_rescale = []
        self.images = []
        self.mapping_filename = None
        self.mfilters = []
        self.mfilters_z = []
        self.n_channels = len(psf_objects)
        self.psf_objects = psf_objects
        self.variances = []
        self.vfilters = []
        self.xt = []
        self.yt = []
        self.z_values = []
            
        # Assert that all the PSF objects are the same size.
        for i in range(1, len(self.psf_objects)):
            assert(self.psf_objects[0].getSize() == self.psf_objects[i].getSize())

        # Update margin based on the psf object size. Note the assumption
        # that all of the psf objects are the same size, or at least smaller
        # than the psf objects for plane 0.
        #
        self.margin = self.psf_objects[0].getMargin()

        # Load the plane to plane mapping data & create affine transform objects.
        mappings = {}
        if parameters.hasAttr("mapping"):
            if os.path.exists(parameters.getAttr("mapping")):
                self.mapping_filename = parameters.getAttr("mapping")
                with open(parameters.getAttr("mapping"), 'rb') as fp:
                    mappings = pickle.load(fp)
            else:
                raise Exception("Mapping file " + parameters.getAttr("mapping") + " does not exist.")

        # Use self.margin - 1, because we added 1 to the x,y coordinates when we saved them.
        for i in range(self.n_channels-1):
            self.xt.append(mpUtil.marginCorrect(mappings["0_" + str(i+1) + "_x"], self.margin - 1))
            self.yt.append(mpUtil.marginCorrect(mappings["0_" + str(i+1) + "_y"], self.margin - 1))
            self.atrans.append(affineTransformC.AffineTransform(xt = self.xt[i],
                                                                yt = self.yt[i]))

        # Note the assumption that the splines for each plane all use
        # the same z scale / have the same z range.
        #
        # 'z_value' is in units of microns.
        #
        self.mfilters_z = parameters.getAttr("z_value", [0.0])
        for zval in self.mfilters_z:
            assert self.psf_objects[0].isValidZ(zval)
            self.z_values.append(self.psf_objects[0].getScaledZ(zval))

        # Configure maxima finder.
        #
        self.mfinder = iaUtilsC.MaximaFinder(margin = self.margin,
                                             radius = self.find_max_radius,
                                             threshold = self.threshold,
                                             z_values = self.z_values)
                                             
        # Load pre-specified peak locations, if any.
        #
        if parameters.hasAttr("peak_locations"):
            [self.peak_locations, self.peak_locations_type] = fitting.getPeakLocations(parameters.getAttr("peak_locations"),
                                                                                       self.margin,
                                                                                       parameters.getAttr("pixel_size"),
                                                                                       self.sigma)

            # Set initial z value (for text files).
            if is_text:
                self.peak_locations["z"][:] = self.z_value[0]

            # Convert z value to PSF FFT units (Insight3 localization files).
            else:
                self.peak_locations["z"] = self.psf_objects[0].getScaledZ(self.peak_locations["z"])

        #
        # Note: A lot of additional initialization occurs in the setVariances() method
        #       because this is the first place where we actually know the image size.
        #
        
    def cleanUp(self):
        super(MPPeakFinder, self).cleanUp()

        # Clean up transforms.
        for at in self.atrans:
            if at is not None:
                at.cleanup()

        # Clean up foreground filters.
        for i in range(len(self.mfilters)):
            for j in range(len(self.mfilters[i])):
                self.mfilters[i][j].cleanup()
    
    def newImage(self, new_images):
        """
        This is called once at the start of the analysis of a new set of images.
        
        new_images - A list of 2D numpy arrays.
        """
        assert (len(new_images) == self.n_channels)
        
        # Initialize new peak minimum threshold.
        #
        if(self.iterations>4):
            self.cur_threshold = self.threshold + 4.0
        else:
            self.cur_threshold = self.threshold + float(self.iterations)

        # Reset maxima finder.
        self.mfinder.resetTaken()

        # Save reference to images.
        #
        self.images = new_images

        # For checking that we're doing the transform correctly and / or have
        # the correct transform.
        if self.check_mode:
            at_images = []
            for i in range(self.n_channels):
                if self.atrans[i] is None:
                    at_images.append(new_images[i].copy())
                else:
                    at_images.append(self.atrans[i].transform(new_images[i]))

            with tifffile.TiffWriter("transform.tif") as tf:
                for at_image in at_images:
                    tf.save(numpy.transpose(at_image.astype(numpy.float32)))

    def peakFinder(self, fit_peaks_images):
        """
        This method does the actual peak finding. 

        We are assuming that Poisson statistics applies, so if we have
        a good estimate of the background and foreground then we can
        predict the significance of a foreground value with reasonable
        accuracy.
        """
        #
        # Calculate (estimated) background variance for each plane.
        #
        # The estimated background and variance should both be > 0.0,
        # or there is going to be trouble.
        #
        bg_variances = []

        # Save fit images for debugging purposes.
        if self.check_mode:
            with tifffile.TiffWriter("fit_images.tif") as tf:
                for fi in fit_images:
                    tf.save(numpy.transpose(fi.astype(numpy.float32)))

        # Iterate over z values.
        for i in range(len(self.vfilters)):
            bg_variance = numpy.zeros(fit_peaks_images[0].shape)

            # Iterate over channels / planes.
            for j in range(len(self.vfilters[i])):

                # Convolve fit image + background with the appropriate variance filter.
                #
                # I believe that this is correct, the variance of the weighted average
                # of independent processes is calculated using the square of the weights.
                #
                conv_var = self.vfilters[i][j].convolve(fit_peaks_images[j] + self.backgrounds[j])

                # Transform variance to the channel 0 frame.
                if self.atrans[j] is None:
                    bg_variance += conv_var
                else:
                    bg_variance += self.atrans[j].transform(conv_var)

            # Camera variances are already convolved and transformed so we just add them on.
            bg_variances.append(bg_variance + self.variances[i])

        # Check for problematic values.
        if self.check_mode:
            for bg in bg_variances:
                mask = (bg <= 0.0)
                if (numpy.sum(mask) > 0):
                    print("Warning! 0.0 / negative values detected in background variance.")
        
        # Save results if needed for debugging purposes.
        if self.check_mode:
            with tifffile.TiffWriter("variances.tif") as tf:
                for bg in bg_variances:
                    tf.save(numpy.transpose(bg.astype(numpy.float32)))
                    
        #
        # Calculate foreground for each z plane.
        #
        fg_averages = []  # This is the average foreground across all the planes for each z value.
        foregrounds = []  # This is the foreground for each plane and z value.

        # Iterate over z values.
        for i in range(len(self.mfilters)):
            foreground = numpy.zeros(fit_peaks_images[0].shape)
            foregrounds.append([])

            # Iterate over channels / planes.
            for j in range(len(self.mfilters[i])):

                # Convolve image / background with the appropriate PSF.
                conv_fg = self.mfilters[i][j].convolve(self.images[j] - fit_peaks_images[j] - self.backgrounds[j])

                # Store convolved image in foregrounds.
                foregrounds[i].append(conv_fg)

                # Transform image to the channel 0 frame.
                if self.atrans[j] is None:
                    foreground += conv_fg
                else:
                    foreground += self.atrans[j].transform(conv_fg)

            fg_averages.append(foreground)

        # Normalize average foreground by background standard deviation.
        #
        fg_bg_ratios = []
        for i in range(len(fg_averages)):
            fg_bg_ratios.append(fg_averages[i]/numpy.sqrt(bg_variances[i]))

        # Save results if needed for debugging purposes.
        #
        if self.check_mode:
            with tifffile.TiffWriter("foregrounds.tif") as tf:
                for fg in fg_averages:
                    tf.save(numpy.transpose(fg.astype(numpy.float32)))

            with tifffile.TiffWriter("fg_bg_ratio.tif") as tf:
                for fg_bg_ratio in fg_bg_ratios:
                    tf.save(numpy.transpose(fg_bg_ratio.astype(numpy.float32)))

        # Apply AOI mask to the images.
        #
        masked_images = []
        for fg_bg_ratio in fg_bg_ratios:
            masked_images.append(fg_bg_ratio * self.peak_mask)

        # Identify local maxima in the masked ratio images stack.
        #
        [x, y, z] = self.mfinder.findMaxima(masked_images)
        return {"x" : x, "y" : y, "z" : z}

    def setVariances(self, variances):

        # Make sure that the number of (sCMOS) variance arrays
        # matches the number of image planes.
        #
        assert(len(variances) == self.n_channels)

        # Pad variances to correct size.
        #
        temp = []
        for variance in variances:
            temp.append(fitting.padArray(variance, self.margin))
        variances = temp
        
        # We initialize the following here because at __init__ we
        # don't know how big the images are.
        #
        # Note the assumption that every frame in all the movies
        # is the same size.
        #
        
        # Create mask to limit peak finding to a user defined sub-region of the image.
        self.peak_mask = numpy.ones(variances[0].shape)
        if self.parameters.hasAttr("x_start"):
            self.peak_mask[0:self.parameters.getAttr("x_start")+self.margin,:] = 0.0
        if self.parameters.hasAttr("x_stop"):
            self.peak_mask[self.parameters.getAttr("x_stop")+self.margin:-1,:] = 0.0
        if self.parameters.hasAttr("y_start"):
            self.peak_mask[:,0:self.parameters.getAttr("y_start")+self.margin] = 0.0
        if self.parameters.hasAttr("y_stop"):
            self.peak_mask[:,self.parameters.getAttr("y_stop")+self.margin:-1] = 0.0

        #
        # Now that we have the MpUtil object we can split the input peak
        # locations to create a list for each channel.
        #
        if self.peak_locations is not None:
            self.peak_locations = self.mpu.splitPeaks(self.peak_locations)

        #
        # Create "foreground" and "variance" filters, as well as the
        # height rescaling array.
        #
        # These are stored in a list indexed by z value, then by
        # channel / plane. So self.mfilters[1][2] is the filter
        # for z value 1, plane 2.
        #
        for i, mfilter_z in enumerate(self.mfilters_z):
            self.height_rescale.append([])
            self.mfilters.append([])
            self.vfilters.append([])

            for j, psf_object in enumerate(self.psf_objects):
                psf = psf_object.getPSF(mfilter_z,
                                        shape = variances[0].shape,
                                        normalize = False)

                #
                # We are assuming that the psf has no negative values,
                # or if it does that they are very small.
                #
                psf_norm = psf/numpy.sum(psf)
                self.mfilters[i].append(matchedFilterC.MatchedFilter(psf_norm))
                self.vfilters[i].append(matchedFilterC.MatchedFilter(psf_norm * psf_norm))

                #
                # This is used to convert the height measured in the
                # convolved image to the correct height in the original
                # image, as this is the height unit that is used in
                # fitting.
                #
                # If you convolved a localization with itself the final
                # height would be sum(loc * loc). Here we are convolving
                # the localizations with a unit sum PSF, so the following
                # should give us the correct initial height under the
                # assumption that the shape of the localization is
                # pretty close to the shape of the PSF.
                #
                self.height_rescale[i].append(1.0/numpy.sum(psf * psf_norm))

                # Save a pictures of the PSFs for debugging purposes.
                if self.check_mode:
                    print("psf max", numpy.max(psf))
                    filename = "psf_z{0:.3f}_c{1:d}.tif".format(mfilter_z, j)
                    tifffile.imsave(filename, numpy.transpose(psf.astype(numpy.float32)))

        # Create matched filter for background.
        bg_psf = fitting.gaussianPSF(variances[0].shape, self.parameters.getAttr("background_sigma"))
        self.bg_filter = matchedFilterC.MatchedFilter(bg_psf)

        #
        # Process variance arrays now as they don't change from frame
        # to frame.
        #
        # This initializes the self.variances array with a list
        # of lists with the same organization as foreground and
        # psf / variance filters.
        #
        # Use variance filter. I now think this is correct as this is
        # also what we are doing with the image background term. In
        # the case of the image background we are estimating the
        # variance under the assumption that it is Poisson so the
        # mean of the background is the variance. With the cameras
        # we know what the variance is because we measured it. Now
        # we need to weight it properly given the PSF filter that
        # we are applying to the foreground.
        #
            
        # Iterate over z values.
        for i in range(len(self.mfilters)):
            variance = numpy.zeros(variances[0].shape)

            # Iterate over channels / planes.
            for j in range(len(self.mfilters[i])):

                # Convolve variance with the appropriate variance filter.
                conv_var = self.vfilters[i][j].convolve(variances[j])

                # Transform variance to the channel 0 frame.
                if self.atrans[j] is None:
                    variance += conv_var
                else:
                    variance += self.atrans[j].transform(conv_var)

            self.variances.append(variance)

        # Save results if needed for debugging purposes.
        if self.check_mode:
            with tifffile.TiffWriter("camera_variances.tif") as tf:
                for var in self.variances:
                    tf.save(numpy.transpose(var.astype(numpy.float32)))

        # Return padded variances
        return variances

    def subtractBackground(self, images, fit_peaks_images, bg_estimates):
        """
        Estimate the background for the images.
        """
        assert(len(images) == self.n_channels)

        if bg_estimates[0] is not None:
            self.backgrounds = bg_estimates

        else:
            self.backgrounds = []
            for i in range(len(images)):
                self.backgrounds.append(self.backgroundEstimator(images[i] - fit_peaks_images[i]))

        # Save results if requested.
        if self.check_mode:
            with tifffile.TiffWriter("bg_estimate.tif") as tf:
                for bg in self.backgrounds:
                    tf.save(numpy.transpose(bg.astype(numpy.float32)))

        return self.backgrounds


class MPPeakFitter(fitting.PeakFitterArbitraryPSF):
    """
    Multi-plane peak fitting.
    """
    def getPeakProperty(self, pname, channel = 0):
        return self.mfitter.getPeakProperty(pname, channel = channel)


class MPFinderFitter(fitting.PeakFinderFitter):
    """
    Multi-plane peak finding and fitting.
    """
    def getPeakProperties(self):
        """
        Create a list of dictionaries with the requested properties.
        """
        peaks = []
        for i in range(self.peak_finder.n_channels):
            temp = {"channel" : i}
            for pname in self.properties:
                temp[pname] = self.peak_fitter.getPeakProperty(pname, channel = i)

                # x,y,z corrections.
                if (pname == "x"):
                    temp[pname] -= float(self.peak_finder.margin)
                    
                elif (pname == "y"):
                    temp[pname] -= float(self.peak_finder.margin)
                
                elif (pname == "z"):
                    temp[pname] = self.peak_fitter.rescaleZ(temp[pname])

            peaks.append(temp)

        return peaks
    
    def loadBackgroundEstimate(self, movie_reader):
        bg_estimates = []
        for i in range(self.peak_finder.n_channels):

            # Load the background of a single channel / plane.
            bg = movie_reader.getBackground(i)

            if bg is None:
                bg_estimates.append(bg)
                continue

            # Add edge padding.
            bg = fitting.padArray(bg, self.peak_finder.margin)

            bg_estimates.append(bg)

        return bg_estimates

    def loadImage(self, movie_reader):
        fit_peaks_images = []
        images = []
        for i in range(self.peak_finder.n_channels):

            # Load the image of a single channel / plane.
            image = movie_reader.getFrame(i)

            # Add edge padding.
            image = fitting.padArray(image, self.peak_finder.margin)

            # Add to lists.
            images.append(image)            
            fit_peaks_images.append(numpy.zeros(image.shape))

        return [images, fit_peaks_images]

    
def initFitter(margin, parameters, psf_objects, variances):
    """
    Create and return a mpFitC.MPSplineFit object.
    """
    assert(len(psf_objects) == len(variances))
    #
    # FIXME: Not sure having two z ranges, one from the spline
    #        and one in the parameters is a good idea.
    #
    # Create the fitter object which will do the actual fitting. Unless specified
    # the fit for each channel is forced to have the same height.
    #
    if isinstance(psf_objects[0], psfFn.PSFFn):
        mfitter = mpFitC.MPPSFFnFit(independent_heights = parameters.getAttr("independent_heights", 0),
                                    psf_objects = psf_objects,
                                    scmos_cals = variances)

    elif isinstance(psf_objects[0], pupilFn.PupilFunction):
        mfitter = mpFitC.MPPupilFnFit(independent_heights = parameters.getAttr("independent_heights", 0),
                                      psf_objects = psf_objects,
                                      scmos_cals = variances)

    elif isinstance(psf_objects[0], splineToPSF.SplineToPSF3D):
        mfitter = mpFitC.MPSplineFit(independent_heights = parameters.getAttr("independent_heights", 0),
                                     psf_objects = psf_objects,
                                     scmos_cals = variances)

    # Pass variances to the fitting object.
    #
    for i in range(len(variances)):
        mfitter.setVariance(variances[i], i)

    # Load mappings.
    #
    if parameters.hasAttr("mapping"):
        if os.path.exists(parameters.getAttr("mapping")):
            mapping_filename = parameters.getAttr("mapping")
        else:
            raise Exception("Mapping file", parameters.getAttr("mapping"), "does not exist.")

        # Use margin - 1, because we added 1 to the x,y coordinates when we saved them.
        #
        mfitter.setMapping(*mpUtil.loadMappings(mapping_filename, margin - 1))
        
    # Load channel Cramer-Rao weights if available.
    #
    weights = None
    if parameters.hasAttr("weights"):
        with open(parameters.getAttr("weights"), 'rb') as fp:
            weights = pickle.load(fp)
    mfitter.setWeights(weights)
    
    return mfitter


def initPSFObjects(parameters):
    """
    Create and return the PSF objects (spline, pupil function or psf FFT).
    """
    psf_objects = []

    # Try PSF FFT.
    #
    if (len(mpUtil.getPSFFFTAttrs(parameters)) > 0):
        
        # Create PSF FFT PSF objects.
        for psf_fft_attr in mpUtil.getPSFFFTAttrs(parameters):
            psf_objects.append(psfFn.PSFFn(psf_filename = parameters.getAttr(psf_fft_attr)))

        # All the PSF FFT objects have to have the same Z range.
        for i in range(1, len(psf_objects)):
            assert (psf_objects[0].getZMin() == psf_objects[i].getZMin())
            assert (psf_objects[0].getZMax() == psf_objects[i].getZMax())
            
    # Try pupil functions.
    #
    elif (len(mpUtil.getPupilFnAttrs(parameters)) > 0):

        # Get fitting Z range (this is in microns).
        [min_z, max_z] = parameters.getZRange()
        
        # Create pupil function PSF objects.
        for pupil_fn_attr in mpUtil.getPupilFnAttrs(parameters):
            psf_objects.append(pupilFn.PupilFunction(pf_filename = parameters.getAttr(pupil_fn_attr),
                                                     zmin = min_z * 1000.0,
                                                     zmax = max_z * 1000.0))

    # Try splines.
    #
    elif (len(mpUtil.getSplineAttrs(parameters)) > 0):

        # Create Spline PSF objects.
        for spline_attr in mpUtil.getSplineAttrs(parameters):
            psf_objects.append(splineToPSF.SplineToPSF3D(spline_file = parameters.getAttr(spline_attr)))

        # All the splines have to have the same Z range.
        for i in range(1, len(psf_objects)):
            assert (psf_objects[0].getZMin() == psf_objects[i].getZMin())
            assert (psf_objects[0].getZMax() == psf_objects[i].getZMax())

    else:
        raise Exception("No PSF objects found.")

    return psf_objects


def initFindAndFit(parameters):
    """
    Create and return a MPFinderFitter object.
    """
    # Create PSF objects.
    psf_objects = initPSFObjects(parameters)
    
    # Create peak finder.
    finder = MPPeakFinder(parameters = parameters,
                          psf_objects = psf_objects)

    # Load sCMOS calibration data.
    #
    # Note: Gain is expected to be in units of ADU per photo-electron.
    #
    variances = []
    for calib_name in mpUtil.getCalibrationAttrs(parameters):
        [offset, variance, gain] = numpy.load(parameters.getAttr(calib_name))
        variances.append(variance/(gain*gain))

    # Set variance in the peak finder. This method also pads the variance
    # to the correct size and performs additional initializations.
    #
    variances = finder.setVariances(variances)

    # Create mpFitC.MPFit object.
    #
    mfitter = initFitter(finder.margin,
                         parameters,
                         psf_objects,
                         variances)

    # Create peak fitter.
    #
    fitter = MPPeakFitter(mfitter = mfitter,
                          parameters = parameters)

    # Specify which properties we want (for each channel) from the
    # analysis. Note that some of these may be duplicates of each
    # other, for example if the heights are not independent.
    #
    properties = ["background", "error", "height", "iterations", "sum", "x", "y", "z"]    

    return MPFinderFitter(peak_finder = finder,
                          peak_fitter = fitter,
                          properties = properties)
