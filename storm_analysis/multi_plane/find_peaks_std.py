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
import storm_analysis.multi_plane.mp_utilities_c as mpUtilC

import storm_analysis.sa_library.affine_transform_c as affineTransformC
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC

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
        self.check_mode = False
        self.height_rescale = []
        self.images = []
        self.mapping_filename = None
        self.mfilters = []
        self.mfilters_z = []
        self.mpu = None
        self.n_channels = len(psf_objects)
        self.psf_objects = psf_objects
        self.variances = []
        self.vfilters = []
        self.xt = []
        self.yt = []
        self.z_values = []

        # Print warning about check mode
        if self.check_mode:
            print("Warning: Running in check mode.")
            
        # Assert that all the PSF objects are the same size.
        for i in range(1, len(self.psf_objects)):
            assert(self.psf_objects[0].getSize() == self.psf_objects[i].getSize())

        #
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
            self.xt.append(mpUtilC.marginCorrect(mappings["0_" + str(i+1) + "_x"], self.margin - 1))
            self.yt.append(mpUtilC.marginCorrect(mappings["0_" + str(i+1) + "_y"], self.margin - 1))
            self.atrans.append(affineTransformC.AffineTransform(xt = self.xt[i],
                                                                yt = self.yt[i]))

        #
        # Note the assumption that the splines for each plane all use
        # the same z scale / have the same z range.
        #
        # 'z_value' is in units of microns.
        #
        self.mfilters_z = parameters.getAttr("z_value", [0.0])
        for zval in self.mfilters_z:
            self.z_values.append(self.psf_objects[0].getScaledZ(zval))

        #
        # Load pre-specified peak locations, if any.
        #
        if parameters.hasAttr("peak_locations"):
            [self.peak_locations, is_text] = fitting.getPeakLocations(parameters.getAttr("peak_locations"),
                                                                      self.margin,
                                                                      parameters.getAttr("pixel_size"),
                                                                      self.sigma)

            zc_index = utilC.getZCenterIndex()
            # Set initial z value (for text files).
            if is_text:
                self.peak_locations[:,zc_index] = self.z_value[0]

            # Convert z value to PSF FFT units (Insight3 localization files).
            else:
                for i in range(self.peak_locations.shape[0]):
                    self.peak_locations[i,zc_index] = self.psf_objects[0].getScaledZ(self.peak_locations[i,zc_index])

        #
        # Note: A lot of additional initialization occurs in the setVariances() method
        #       because this is the first place where we actually know the image size.
        #
        
    def cleanUp(self):
        self.mpu.cleanup()
        for at in self.atrans:
            if at is not None:
                at.cleanup()

    def getMPU(self):
        """
        Return the MPU object (this object handles manipulation of the peak lists).
        """
        return self.mpu

    def mergeNewPeaks(self, peaks, new_peaks):
        return self.mpu.mergeNewPeaks(peaks, new_peaks)
    
    def newImages(self, new_images):
        """
        This is called once at the start of the analysis of a new set of images.
        
        new_images - A list of 2D numpy arrays.
        """
        #
        # Initialize new peak minimum threshold.
        #
        if(self.iterations>4):
            self.cur_threshold = self.threshold + 4.0
        else:
            self.cur_threshold = self.threshold + float(self.iterations)

        #
        # Reset taken arrays.
        #
        self.taken = []
        for i in range(len(self.mfilters_z)):
            self.taken.append(numpy.zeros(new_images[0].shape, dtype=numpy.int32))

        #
        # Save references to images & create empty list for the background estimates.
        #
        self.images = new_images
        self.backgrounds = []
        for i in range(len(self.images)):
            self.backgrounds.append(None)

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

    def peakFinder(self, fit_images):
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
            bg_variance = numpy.zeros(fit_images[0].shape)

            # Iterate over channels / planes.
            for j in range(len(self.vfilters[i])):

                # Convolve fit image + background with the appropriate variance filter.
                #
                # I believe that this is correct, the variance of the weighted average
                # of independent processes is calculated using the square of the weights.
                #
                conv_var = self.vfilters[i][j].convolve(fit_images[j] + self.backgrounds[j])

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
            foreground = numpy.zeros(fit_images[0].shape)
            foregrounds.append([])

            # Iterate over channels / planes.
            for j in range(len(self.mfilters[i])):

                # Convolve image / background with the appropriate PSF.
                conv_fg = self.mfilters[i][j].convolve(self.images[j] - fit_images[j] - self.backgrounds[j])

                # Store convolved image in foregrounds.
                foregrounds[i].append(conv_fg)

                # Transform image to the channel 0 frame.
                if self.atrans[j] is None:
                    foreground += conv_fg
                else:
                    foreground += self.atrans[j].transform(conv_fg)

            fg_averages.append(foreground)

        # Normalize average foreground by background standard deviation.
        fg_bg_ratios = []
        for i in range(len(fg_averages)):
            fg_bg_ratios.append(fg_averages[i]/numpy.sqrt(bg_variances[i]))

        # Save results if needed for debugging purposes.
        if self.check_mode:
            with tifffile.TiffWriter("foregrounds.tif") as tf:
                for fg in fg_averages:
                    tf.save(numpy.transpose(fg.astype(numpy.float32)))

            with tifffile.TiffWriter("fg_bg_ratio.tif") as tf:
                for fg_bg_ratio in fg_bg_ratios:
                    tf.save(numpy.transpose(fg_bg_ratio.astype(numpy.float32)))

        #
        # At each z value, find peaks in foreground image normalized
        # by the background standard deviation.
        #
        all_new_peaks = None
        zero_array = numpy.zeros(fg_bg_ratios[0].shape)
        for i in range(len(self.mfilters)):

            #
            # Mask the image so that peaks are only found in the AOI. Ideally the
            # this mask should probably be adjusted to limit analysis to only
            # the regions of the image where there is data from every channel / plane.
            #
            masked_image = fg_bg_ratios[i] * self.peak_mask
        
            # Identify local maxima in the masked image.
            [new_peaks, taken] = utilC.findLocalMaxima(masked_image,
                                                       self.taken[i],
                                                       self.cur_threshold,
                                                       self.find_max_radius,
                                                       self.margin)

            #
            # Initialize peaks with normalized height value. We'll split these
            # later into peaks for each plane, and at that point the height,
            # background and z values will be corrected.
            #
            # Note: Sigma is irrelevant for fitting, but it needs to be some non-zero number.
            #
            new_peaks = utilC.initializePeaks(new_peaks,    # The new peaks.
                                              masked_image, # Use SNR as height, corrected later for fitting.
                                              zero_array,   # Zero for now, corrected later for fitting.
                                              self.sigma,   # The starting sigma value.
                                              i)            # Index of the z-plane, the actual z value is added later.

            # Add to all peaks accumulator.
            if all_new_peaks is None:
                all_new_peaks = new_peaks
            else:
                all_new_peaks = numpy.append(all_new_peaks, new_peaks, axis = 0)

        #
        # If there are multiple peaks with similar x,y but in different
        # planes, use the one with the highest normalized value.
        #
        # FIXME: If the planes are far enough apart in z we should allow
        #        peaks with a similar x,y.
        #
        if (len(self.mfilters) > 1):
            all_new_peaks = utilC.removeClosePeaks(all_new_peaks,                                               
                                                   self.find_max_radius,
                                                   self.find_max_radius)

        #
        # Split into a peak/localization for each image plane.
        #
        # Note that the peaks array is expected to have all the peaks
        # for the first plane first, then all the peaks for the second
        # plane, etc.. With the same number of peaks per plane.
        #
        # This is how you would access the same peak in different channels:
        #
        # ch0) all_new_peaks[0 * n_peaks + peak_number]
        # ch1) all_new_peaks[1 * n_peaks + peak_number]
        # etc..
        #
        all_new_peaks = self.mpu.splitPeaks(all_new_peaks)

        #
        # Remove peaks with members in one or more channels that are
        # outside of the image.
        #
        all_new_peaks = self.mpu.filterPeaks(all_new_peaks, self.mpu.badPeakMask(all_new_peaks))

        # Initialize background values.
        mpUtilC.initializeBackground(all_new_peaks, self.backgrounds)

        # Need to do this before z initialization as we are using the
        # z value to index into the foregrounds array.
        mpUtilC.initializeHeight(all_new_peaks, foregrounds, self.height_rescale)

        # Replace z index with the z value used as the initial guess
        # for fitting.
        mpUtilC.initializeZ(all_new_peaks, self.z_values)
        
        if False:
            pp = 3
            if (all_new_peaks.shape[0] > pp):
                for i in range(pp):
                    print("Peak",i)
                    self.mpu.prettyPrintPeak(all_new_peaks, i)
                    print("")
        
        return all_new_peaks

    def setVariances(self, variances):

        #
        # Make sure that the number of (sCMOS) variance arrays
        # matches the number of image planes.
        #
        assert(len(variances) == self.n_channels)

        #
        # Pad variances to correct size.
        #
        temp = []
        for variance in variances:
            temp.append(fitting.padArray(variance, self.margin))
        variances = temp
        
        #
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
        # Create mpUtilC.MpUtil object that is used to do a lot of the
        # peak list manipulations.
        #
        self.mpu = mpUtilC.MpUtil(radius = self.new_peak_radius,
                                  neighborhood = self.neighborhood,
                                  im_size_x = variances[0].shape[1],
                                  im_size_y = variances[0].shape[0],
                                  n_channels = self.n_channels,
                                  n_zplanes = len(self.z_values),
                                  margin = self.margin)

        #
        # Load mappings file again so that we can set the transforms for
        # the MpUtil object.
        #
        # Use self.margin - 1, because we added 1 to the x,y coordinates
        # when we saved them, see sa_library.i3dtype.createFromMultiFit().
        #
        [xt, yt] = mpUtilC.loadMappings(self.mapping_filename, self.margin - 1)[:2]
        self.mpu.setTransforms(xt, yt)

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
                if False:
                    print("psf max", numpy.max(psf))
                    filename = "psf_z{0:.3f}_c{1:d}.tif".format(mfilter_z, j)
                    tifffile.imsave(filename, psf.astype(numpy.float32))

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

    def subtractBackground(self, image, bg_estimate, index):
        """
        Estimate the background for the image from a particular
        plane (specified by index).

        Note: image is the residual image after the found / fit
              localizations have been subtracted out.
        """
        if bg_estimate is not None:
            self.backgrounds[index] = bg_estimate

        else:
            self.backgrounds[index] = self.backgroundEstimator(image)

        # Save results if needed for debugging purposes.
        if self.check_mode and (index == (self.n_channels - 1)):
            with tifffile.TiffWriter("bg_estimate.tif") as tf:
                for bg in self.backgrounds:
                    tf.save(numpy.transpose(bg.astype(numpy.float32)))


class MPPeakFitter(fitting.PeakFitter):
    """
    Multi-plane peak fitting.
    """
    def __init__(self, mpu = None, n_channels = None, **kwds):
        super(MPPeakFitter, self).__init__(**kwds)
        self.images = None
        self.mpu = mpu
        self.n_channels = n_channels

    def cleanUp(self):
        super(MPPeakFitter, self).cleanUp()
        self.mpu.cleanup()

    def fitPeaks(self, peaks):
        
        # Fit to update peak locations.
        [fit_peaks, fit_peaks_images] = self.peakFitter(peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks, 0.0)

        # Remove peaks that are too close to each.
        fit_peaks = self.mpu.removeClosePeaks(fit_peaks)

        # Update fits for remaining peaks.
        #
        # FIXME: Check if it makes sense do this if haven't removed any
        #        peaks. Do we need the extra iterations?
        #
        [fit_peaks, fit_peaks_images] = self.peakFitter(fit_peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks, 0.0)

        # Save fit images for debugging.
        if False:
            with tifffile.TiffWriter("fit_images.tif") as tf:
                for fp_im in fit_peaks_images:
                    tf.save(numpy.transpose(fp_im.astype(numpy.float32)))
        
        return [fit_peaks, fit_peaks_images]

    def newImages(self, new_images):
        for i, image in enumerate(new_images):
            self.mfitter.newImage(image, i)

    def peakFitter(self, peaks):
        """
        This method does the actual peak fitting.
        """
        fit_peaks = self.mfitter.doFit(peaks)
        fit_peaks_images = []
        for i in range(self.n_channels):
            fit_peaks_images.append(self.mfitter.getFitImage(i))
        return [fit_peaks, fit_peaks_images]

    def rescaleZ(self, peaks):
        """
        Convert from spline z units to real z units.
        """
        return self.mfitter.rescaleZ(peaks)


class MPFinderFitter(fitting.PeakFinderFitter):
    """
    Multi-plane peak finding and fitting.
    """
    def __init__(self, n_planes = None, **kwds):
        super(MPFinderFitter, self).__init__(**kwds)
        
        self.n_planes = n_planes

    def analyzeImage(self, movie_reader, save_residual = False, verbose = False):
        """
        Analyze an "image" and return a list of the found localizations.

        Internally the list of localizations is a multiple of the number of planes.
        """
        #
        # Load and scale the images.
        #
        [images, fit_peaks_images] = self.loadImages(movie_reader)

        #
        # Load and scale the background estimates.
        #
        bg_estimates = self.loadBackgroundEstimates(movie_reader)
        
        self.peak_finder.newImages(images)
        self.peak_fitter.newImages(images)

        peaks = False
        for i in range(self.peak_finder.iterations):
            if verbose:
                print(" iteration", i)

            # Update background estimate.
            for j in range(self.n_planes):
                self.peak_finder.subtractBackground(images[j] - fit_peaks_images[j], bg_estimates[j], j)

            # Find new peaks.
            [found_new_peaks, peaks] = self.peak_finder.findPeaks(fit_peaks_images, peaks)
            
            # Fit new peaks.
            if isinstance(peaks, numpy.ndarray):
                if verbose:
                    print("  found", peaks.shape[0]/self.n_planes)
                [peaks, fit_peaks_images] = self.peak_fitter.fitPeaks(peaks)
                if verbose:
                    print("  fit", peaks.shape[0]/self.n_planes)

            if not found_new_peaks:
                break

        if isinstance(peaks, numpy.ndarray):
            peaks[:,utilC.getXCenterIndex()] -= float(self.margin)
            peaks[:,utilC.getYCenterIndex()] -= float(self.margin)

        #
        # sa_utilities.std_analysis doesn't do anything with the second
        # argument, historically the residual, so just return None.
        #
        return [peaks, None]

    def cleanUp(self):
        self.peak_finder.cleanUp()
        self.peak_fitter.cleanUp()

    def getConvergedPeaks(self, peaks):
        converged_peaks = super().getConvergedPeaks(peaks)
        return self.peak_fitter.rescaleZ(converged_peaks)

    def loadBackgroundEstimates(self, movie_reader):
        bg_estimates = []
        for i in range(self.n_planes):

            # Load the background of a single channel / plane.
            bg = movie_reader.getBackground(i)

            if bg is None:
                bg_estimates.append(bg)
                continue

            # Add edge padding.
            bg = fitting.padArray(bg, self.margin)

            bg_estimates.append(bg)

        return bg_estimates
        
    def loadImages(self, movie_reader):
        fit_peaks_images = []
        images = []
        for i in range(self.n_planes):

            # Load the image of a single channel / plane.
            image = movie_reader.getFrame(i)

            # Add edge padding.
            image = fitting.padArray(image, self.margin)

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
    if isinstance(psf_objects[0], pupilFn.PupilFunction):
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
        mfitter.setMapping(*mpUtilC.loadMappings(mapping_filename, margin - 1))
        
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

    # Try pupil functions.
    #
    if (len(mpUtilC.getPupilFnAttrs(parameters)) > 0):

        # Get fitting Z range (this in microns).
        [min_z, max_z] = parameters.getZRange()
        
        # Create pupil function PSF objects.
        for pupil_fn_attr in mpUtilC.getPupilFnAttrs(parameters):
            psf_objects.append(pupilFn.PupilFunction(parameters.getAttr(pupil_fn_attr),
                                                     zmin = min_z * 1000.0,
                                                     zmax = max_z * 1000.0))

    # Try splines.
    #
    elif (len(mpUtilC.getSplineAttrs(parameters)) > 0):

        # Create Spline PSF objects.
        for spline_attr in mpUtilC.getSplineAttrs(parameters):
            psf_objects.append(splineToPSF.SplineToPSF3D(parameters.getAttr(spline_attr)))

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

    # Get margin size from finder.
    margin = finder.margin

    # Load sCMOS calibration data.
    #
    # Note: Gain is expected to be in units of ADU per photo-electron.
    #
    n_planes = 0
    variances = []
    for calib_name in mpUtilC.getCalibrationAttrs(parameters):
        n_planes += 1
        [offset, variance, gain] = numpy.load(parameters.getAttr(calib_name))
        variances.append(variance/(gain*gain))

    # Set variance in the peak finder. This method also pads the variance
    # to the correct size and performs additional initializations.
    #
    variances = finder.setVariances(variances)

    # Create mpFitC.MPFit object.
    #
    mfitter = initFitter(margin, parameters, psf_objects, variances)

    # Create peak fitter.
    #
    fitter = MPPeakFitter(mfitter = mfitter,
                          mpu = finder.getMPU(),
                          n_channels = n_planes,
                          parameters = parameters)

    return MPFinderFitter(peak_finder = finder,
                          peak_fitter = fitter,
                          n_planes = n_planes)
