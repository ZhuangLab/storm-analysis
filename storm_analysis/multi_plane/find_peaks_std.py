#!/usr/bin/env python
"""
(Standard) multi-plane peak finder and fitter.

"channel", "frame", "image" and "plane" are used somewhat interchangeably..

Hazen 05/17
"""

import pickle
import numpy
import tifffile

import storm_analysis.multi_plane.mp_fit_c as mpFitC
import storm_analysis.multi_plane.mp_utilities_c as mpUtilC

import storm_analysis.sa_library.affine_transform_c as affineTransformC
import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC

import storm_analysis.sa_utilities.std_analysis as stdAnalysis

import storm_analysis.simulator.draw_gaussians_c as dg

import storm_analysis.spliner.spline_to_psf as splineToPSF


def getAttrs(parameters, pre, post, max_value = 8):
    pnames = []
    for i in range(max_value):
        pname = pre + str(i) + post
        if parameters.hasAttr(pname):
            pnames.append(pname)
    return pnames

def getCalibrationAttrs(parameters):
    return getAttrs(parameters, "channel", "_cal")

def getExtAttrs(parameters):
    return getAttrs(parameters, "channel", "_ext")

def getOffsetAttrs(parameters):
    return getAttrs(parameters, "channel", "_offset")

def getSplineAttrs(parameters):
    return getAttrs(parameters, "spline", "")


class MPMovieReader(stdAnalysis.MovieReader):
    """
    Movie reader specialized for multi-plane data.
    """
    def __init__(self, base_name = None, parameters = None):

        self.backgrounds = []
        self.bg_estimators = []
        self.cur_frame = 0
        self.frames = []
        self.max_frame = 0
        self.offsets = []
        self.parameters = parameters
        self.planes = []

        # Load the movies and offsets for each plane/channel.
        for ext in getExtAttrs(parameters):
            movie_name = base_name + parameters.getAttr(ext)
            print(movie_name)
            self.planes.append(datareader.inferReader(movie_name))

        for offset in getOffsetAttrs(parameters):
            self.offsets.append(parameters.getAttr(offset))

        print("Found data for", len(self.planes), "planes.")

        [self.movie_x, self.movie_y, self.movie_l] = self.planes[0].filmSize()

    def getBackground(self, plane):
        if (len(self.backgrounds) > 0):
            return self.backgrounds[plane]
        else:
            return None

    def getFrame(self, plane):
        return self.frames[plane]

    def hashID(self):
        return self.planes[0].hashID()

    def nextFrame(self):
        if (self.cur_frame < self.max_frame):

            # Update background estimate.
            self.backgrounds = []
            for i, bg_estimator in enumerate(self.bg_estimators):
                self.backgrounds.append(bg_estimator.estimateBG(self.cur_frame + self.offsets[i]))

            # Load planes. Zero / negative value removal is done later in the
            # analysis (after offset subtraction and gain correction).
            self.frames = []
            for i, plane in enumerate(self.planes):
                self.frames.append(plane.loadAFrame(self.cur_frame + self.offsets[i]))

            self.cur_frame += 1
            return True
        else:
            return False
        
    def setup(self, start_frame):

        # Figure out where to start.
        self.cur_frame = start_frame
        if self.parameters.hasAttr("start_frame"):
            if (self.parameters.getAttr("start_frame") >= self.cur_frame):
                if (self.parameters.getAttr("start_frame") < self.movie_l):
                    self.cur_frame = self.parameters.getAttr("start_frame")
        
        # Figure out where to stop.
        self.max_frame = self.movie_l
        if self.parameters.hasAttr("max_frame"):
            if (self.parameters.getAttr("max_frame") > 0):
                if (self.parameters.getAttr("max_frame") < self.movie_l):
                    self.max_frame = self.parameters.getAttr("max_frame")

        # Configure background estimator, if any.
        if (self.parameters.getAttr("static_background_estimate", 0) > 0):
            print("Using static background estimator.")
            s_size = self.parameters.getAttr("static_background_estimate")
            for i in range(len(self.planes)):
                bg_est = static_background.StaticBGEstimator(self.planes[i],
                                                             start_frame = self.cur_frame + self.offsets[i],
                                                             sample_size = s_size)
                self.bg_estimators.append(bg_est)


class MPPeakFinder(fitting.PeakFinder):
    """
    Multi-plane peak finding.

    All locations are relative to plane 0. The number of localizations
    is a multiple of the number of focal planes.

    The expectation is that we are working with an image that has already 
    been corrected for gain and offset.

    This works with affine transformed versions of the images so that they
    can all be overlaid on top of each other.
    """
    def __init__(self, parameters):
        super().__init__(parameters)

        self.atrans = [None]
        self.backgrounds = []
        self.bg_filter = None
        self.bg_filter_sigma = parameters.getAttr("bg_filter_sigma")
        self.height_rescale = []
        self.images = []
        self.initialized = False
        self.mfilters = []
        self.mfilters_z = []
        self.n_channels = 0
        self.s_to_psfs = []
        self.variances = []
        self.vfilters = []
        self.xt = []
        self.yt = []
        self.z_values = []
                  
        # Load the splines.
        for spline_name in getSplineAttrs(parameters):
            self.s_to_psfs.append(splineToPSF.loadSpline(parameters.getAttr(spline_name)))
            self.n_channels += 1

        # Update margin based on the spline size. Note the assumption
        # that all of the splines are the same size, or at least smaller
        # than spline for plane 0.
        old_margin = self.margin
        self.margin = int((self.s_to_psfs[0].getSize() + 1)/4 + 2)

        # Load the plane to plane mapping data & create affine transform objects.
        mappings = {}
        if parameters.hasAttr("mapping"):
            if os.path.exists(parameters.getAttr("mapping")):
                with open(parameters.getAttr("mapping"), 'rb') as fp:
                    mappings = pickle.load(fp)

        for i in range(self.n_channels-1):
            self.xt.append(mpUtilC.marginCorrect(mappings["0_" + str(i+1) + "_x"], self.margin))
            self.yt.append(mpUtilC.marginCorrect(mappings["0_" + str(i+1) + "_y"], self.margin))
            self.atrans.append(affineTransformC.AffineTransform(xt = self.xt[i],
                                                                yt = self.yt[i]))
        
        # Note the assumption that the splines for each plane all use
        # the same z scale / have the same z range.
        self.mfilters_z = parameters.getAttr("z_value", [0.0])
        for zval in self.mfilters_z:
            self.z_values.append(self.s_to_psfs[0].getScaledZ(zval))

        if parameters.hasAttr("peak_locations"):

            # Correct for any difference in the margins.
            self.peak_locations[:,utilC.getXCenterIndex()] += self.margin - old_margin
            self.peak_locations[:,utilC.getYCenterIndex()] += self.margin - old_margin

            # Provide the "correct" starting z value.
            #
            # FIXME: Should also allow the user to specify the peak z location
            #        as any fixed value could be so far off for some of the
            #        localizations that the fitting will fail.
            #
            self.peak_locations[:,utilC.getZCenterIndex()] = self.z_values[0]

            # Split into per channel peaks.
            self.peak_locations = mpUtilC.splitPeaks(self.peak_locations, self.xt, self.yt)

    def backgroundEstimator(self, image):
        #
        # FIXME: May have ringing if the background is not zero at the edge
        #        of the image.
        #
        return self.bg_filter.convolve(image)
        
    def cleanUp(self):
        for at in self.atrans:
            if at is not None:
                at.cleanup()

    def findPeaks(self, fit_images, peaks):
        """
        Finds the peaks in an image & adds to the current list of peaks.
   
        fit_images - An image for each plane containing the current fit peaks
                     without the background term.
        peaks - The current list of peaks.
    
        return - [True/False if new peaks were added to the current list, the new peaks]
        """
        #
        # This is basically identical to the equivalent method in the fitting.PeakFinder
        # super class, but we need to handle the merge differently as we have peaks
        # for each channel.
        #

        # Use pre-specified peak locations if available, e.g. bead calibration.
        if self.peak_locations is not None:
            new_peaks = self.peak_locations
            
        # Otherwise, identify local maxima in the image and initialize fitting parameters.
        else:
            new_peaks = self.peakFinder(fit_images)

        # Update new peak identification threshold (if necessary).
        # Also, while threshold is greater than min_threshold we
        # are automatically not done.
        found_new_peaks = False
        if (self.cur_threshold > self.threshold):
            self.cur_threshold -= self.threshold
            found_new_peaks = True

        # If we did not find any new peaks then we may be done.
        if (new_peaks.shape[0] == 0):
            return [found_new_peaks, peaks]

        # Add new peaks to the current list of peaks if it exists,
        # otherwise these peaks become the current list.
        if isinstance(peaks, numpy.ndarray):
            merged_peaks = mpUtilC.mergeNewPeaks(peaks,
                                                 new_peaks,
                                                 self.new_peak_radius,
                                                 self.neighborhood,
                                                 self.n_channels)
        
            # If none of the new peaks are valid then we may be done.
            if (merged_peaks.shape[0] == peaks.shape[0]):
                return [found_new_peaks, merged_peaks]
            else:
                return [True, merged_peaks]
        else:
            return [True, new_peaks]
    
    def newImages(self, new_images):
        """
        This is called once at the start of the analysis of a new set of images.
        
        new_images - A list of 2D numpy arrays.
        """
        #
        # Initialize new peak minimum threshold.
        #
        if(self.iterations>4):
            self.cur_threshold = 4.0 * self.threshold
        else:
            self.cur_threshold = float(self.iterations) * self.threshold

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
        if True:
            at_images = []
            for i in range(self.n_channels):
                if self.atrans[i] is None:
                    at_images.append(new_images[i].copy())
                else:
                    at_images.append(self.atrans[i].transform(new_images[i]))

            with tifffile.TiffWriter("transform.tif") as tf:
                for at_image in at_images:
                    tf.save(at_image.astype(numpy.float32))

        #
        # We initialize the following here because at __init__ we
        # don't know how big the images are.
        #
        # Note the assumption that every frame in all the movies
        # is the same size.
        #
        if not self.initialized:
            assert(len(new_images) == self.n_channels)
            
            self.initialized = True
        
            # Create mask to limit peak finding to a user defined sub-region of the image.
            self.peak_mask = numpy.ones(new_images[0].shape)
            if self.parameters.hasAttr("x_start"):
                self.peak_mask[0:self.parameters.getAttr("x_start")+self.margin,:] = 0.0
            if self.parameters.hasAttr("x_stop"):
                self.peak_mask[self.parameters.getAttr("x_stop")+self.margin:-1,:] = 0.0
            if self.parameters.hasAttr("y_start"):
                self.peak_mask[:,0:self.parameters.getAttr("y_start")+self.margin] = 0.0
            if self.parameters.hasAttr("y_stop"):
                self.peak_mask[:,self.parameters.getAttr("y_stop")+self.margin:-1] = 0.0

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
                
                for j, s_to_psf in enumerate(self.s_to_psfs):
                    psf = s_to_psf.getPSF(mfilter_z,
                                          shape = new_images[0].shape,
                                          normalize = False)

                    #
                    # We are assuming that the psf has no negative values,
                    # or if it does that they are very small.
                    #
                    psf_norm = psf/numpy.sum(psf)
                    self.mfilters[i].append(matchedFilterC.MatchedFilter(psf_norm))
                    self.vfilters[i].append(matchedFilterC.MatchedFilter(psf_norm * psf_norm))
                    
                    self.height_rescale[i].append(1.0/numpy.sum(psf * psf_norm))

                    # Save a pictures of the PSFs for debugging purposes.
                    if False:
                        print("psf max", numpy.max(psf))
                        filename = "psf_z{0:.3f}_c{1:d}.tif".format(mfilter_z, j)
                        tifffile.imsave(filename, psf.astype(numpy.float32))

            # "background" filter.
            psf = dg.drawGaussiansXY(new_images[0].shape,
                                     numpy.array([0.5*new_images[0].shape[0]]),
                                     numpy.array([0.5*new_images[0].shape[1]]),
                                     sigma = self.bg_filter_sigma)
            psf = psf/numpy.sum(psf)
            self.bg_filter = matchedFilterC.MatchedFilter(psf)

            #
            # Process variance arrays now as they don't change from frame
            # to frame.
            #
            # This converts the original self.variances array into a list
            # of lists with the same organization as foreground and
            # variance filters.
            #
            variances = self.variances
            self.variances = []
            
            # Iterate over z values.
            for i in range(len(self.vfilters)):
                variance = numpy.zeros(variances[0].shape)

                # Iterate over channels / planes.
                for j in range(len(self.vfilters[i])):

                    # Convolve variance with the appropriate variance filter.
                    conv_var = self.vfilters[i][j].convolve(variances[j])

                    # Transform variance to the channel 0 frame.
                    if self.atrans[j] is None:
                        variance += conv_var
                    else:
                        variance += self.atrans[j].transform(conv_var)

                self.variances.append(variance)

    def peakFinder(self, fit_images):
        """
        This method does the actual peak finding.
        """
        #
        # Calculate (estimated) background variance for each plane.
        #
        # The estimated background and variance should both be > 0.0,
        # or there is going to be trouble.
        #
        bg_variances = []

        # Iterate over z values.
        for i in range(len(self.vfilters)):
            bg_variance = numpy.zeros(fit_images[0].shape)

            # Iterate over channels / planes.
            for j in range(len(self.vfilters[i])):

                # Convolve fit image + background with the appropriate variance filter.
                conv_var = self.vfilters[i][j].convolve(fit_images[j] + self.backgrounds[j])

                # Transform variance to the channel 0 frame.
                if self.atrans[j] is None:
                    bg_variance += conv_var
                else:
                    bg_variance += self.atrans[j].transform(conv_var)

            bg_variances.append(bg_variance + self.variances[i])

        # Check for problematic values.
        if True:
            for bg in bg_variances:
                mask = (bg <= 0.0)
                if (numpy.sum(mask) > 0):
                    print("Warning! 0.0 / negative values detected in background variance.")
        
        # Save results if needed for debugging purposes.
        if True:
            with tifffile.TiffWriter("variances.tif") as tf:
                for bg in bg_variances:
                    tf.save(bg.astype(numpy.float32))
                    
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
                conv_fg = self.mfilters[i][j].convolve(self.images[j] - self.backgrounds[j])

                # Store convolved image in foregrounds.
                foregrounds[i].append(conv_fg)

                # Transform image to the channel 0 frame.
                if self.atrans[j] is None:
                    foreground += conv_fg
                else:
                    foreground += self.atrans[j].transform(conv_fg)

            fg_averages.append(foreground * self.peak_mask)

        # Normalize average foreground by background standard deviation.
        fg_bg_ratios = []
        for i in range(len(fg_averages)):
            fg_bg_ratios.append(fg_averages[i]/numpy.sqrt(bg_variances[i]))

        # Save results if needed for debugging purposes.
        if True:
            with tifffile.TiffWriter("foregrounds.tif") as tf:
                for fg in fg_averages:
                    tf.save(fg.astype(numpy.float32))                    

            with tifffile.TiffWriter("fg_bg_ratio.tif") as tf:
                for fg_bg_ratio in fg_bg_ratios:
                    tf.save(fg_bg_ratio.astype(numpy.float32))

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
        all_new_peaks = mpUtilC.splitPeaks(all_new_peaks, self.xt, self.yt)
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
                    mpUtilC.prettyPrintPeak(all_new_peaks, i, self.n_channels)
                    print("")
        
        return all_new_peaks

    def setVariances(self, variances):
        
        # Make sure that the number of (sCMOS) variance arrays
        # matches the number of image planes.
        assert(len(variances) == self.n_channels)

        self.variances = variances

    def subtractBackground(self, image, bg_estimate, index):
        """
        Estimate the background for the image from a particular
        plane (specified by index).

        Note: image is the residual image after the found / fit
              localizations have been subtracted out.

        Note: Unlike PeakFinder.subtractBackground this method
              does not return anything.
        """
        if bg_estimate is not None:
            self.backgrounds[index] = bg_estimate

        else:
            self.backgrounds[index] = self.backgroundEstimator(image)

        # Save results if needed for debugging purposes.
        if True and (index == (self.n_channels - 1)):
            with tifffile.TiffWriter("bg_estimate.tif") as tf:
                for bg in self.backgrounds:
                    tf.save(bg.astype(numpy.float32))


class MPPeakFitter(fitting.PeakFitter):
    """
    Multi-plane peak fitting.
    """
    def __init__(self, parameters):
        super().__init__(parameters)
        self.images = None
        self.mappings = None
        self.margin = 0
        self.n_channels = None
        self.variances = []
        
        # Load the plane to plane mapping.
        self.mappings = {}
        if parameters.hasAttr("mapping"):
            if os.path.exists(parameters.getAttr("mapping")):
                with open(parameters.getAttr("mapping"), 'rb') as fp:
                    self.mappings = pickle.load(fp)
            
        # Load the splines & create the multi-plane spline fitter.
        coeffs = []
        splines = []
        for spline_name in getSplineAttrs(parameters):
            with open(parameters.getAttr(spline_name), 'rb') as fp:
                spline_data = pickle.load(fp)
            self.zmin = spline_data["zmin"]/1000.0
            self.zmax = spline_data["zmax"]/1000.0
            coeffs.append(spline_data["coeff"])
            splines.append(spline_data["spline"])

        self.mfitter = mpFitC.MPSplineFit(splines, coeffs, verbose = False)

        self.n_channels = len(splines)

    def fitPeaks(self, peaks):

        # Fit to update peak locations.
        [fit_peaks, fit_peaks_images] = self.peakFitter(peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks, 0.9 * self.threshold)

        # Save fit images for debugging.
        if True:
            with tifffile.TiffWriter("fit_images.tif") as tf:
                for fp_im in fit_peaks_images:
                    tf.save(fp_im.astype(numpy.float32))

        # Remove peaks that are too close to each other & refit.
        fit_peaks = mpUtilC.removeClosePeaks(fit_peaks, self.sigma, self.neighborhood)
        
        [fit_peaks, fit_peaks_images] = self.peakFitter(fit_peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks, 0.9 * self.threshold)

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
        return self.mfitter.rescaleZ(peaks, self.zmin, self.zmax)

    def setMargin(self, margin):
        self.margin = margin

        # Update margins and pass the fitting object.
        for map_key in self.mappings:
            [ch_from, ch_to, x_or_y] = map_key.split("_")
            mp = mpUtilC.marginCorrect(self.mappings[map_key], self.margin)
            self.mfitter.setMapping(int(ch_from), int(ch_to), mp, (x_or_y == "x"))

    def setVariances(self, variances):
        """
        Set fitter (sCMOS) variance arrays.
        """
        assert(len(variances) == self.n_channels)

        # Pass variances to the fitting object.
        for i in range(self.n_channels):
            self.mfitter.setVariance(variances[i], i)

    def setWeights(self):
        """
        Tells the fitter object to pass the channel and z dependent
        parameter weight values to the C library.
        """
        self.mfitter.setWeights()

            
class MPFinderFitter(fitting.PeakFinderFitter):
    """
    Multi-plane peak finding and fitting.
    """
    def __init__(self, parameters):
        super().__init__(parameters)
        self.gains = []
        self.n_planes = 0
        self.offsets = []
        self.peak_finder = MPPeakFinder(parameters)
        self.peak_fitter = MPPeakFitter(parameters)
        self.variances = []

        # Update margin.
        self.margin = self.peak_finder.margin
        self.peak_fitter.setMargin(self.margin)

        # Load sCMOS calibration data.
        #
        # Note: Gain is expected to be in units of ADU per photo-electron.
        #
        for calib_name in getCalibrationAttrs(parameters):
            [offset, variance, gain] = fitting.loadSCMOSData(parameters.getAttr(calib_name),
                                                             self.margin)
            self.offsets.append(offset)
            self.variances.append(variance/(gain*gain))
            self.gains.append(1.0/gain)

            self.n_planes += 1

        self.peak_finder.setVariances(self.variances)
        self.peak_fitter.setVariances(self.variances)
        self.peak_fitter.setWeights()

    def analyzeImage(self, movie_reader, save_residual = False, verbose = False):
        """
        Analyze an "image" and return a list of the found localizations.

        Internally the list of localizations is a multiple of the number of planes.
        """
        #
        # Load and scale the images.
        #
        fit_peaks_images = []
        images = []
        for i in range(self.n_planes):

            # Load the image of a single channel / plane.
            image = movie_reader.getFrame(i)

            # Add edge padding.
            image = fitting.padArray(image, self.margin)

            # Convert to photo-electrons.
            image = (image - self.offsets[i])*self.gains[i]

            # Remove values < 1.0
            mask = (image < 1.0)
            if (numpy.sum(mask) > 0):
                image[mask] = 1.0

            images.append(image)

            fit_peaks_images.append(numpy.zeros(image.shape))

        #
        # Load and scale the background estimates.
        #
        bg_estimates = []
        for i in range(self.n_planes):

            # Load the background of a single channel / plane.
            bg = movie_reader.getBackground(i)

            if bg is None:
                bg_estimates.append(bg)
                continue

            # Add edge padding.
            bg = fitting.padArray(bg, self.margin)

            # Convert to photo-electrons.
            bg = (bg - self.offsets[i])*self.gains[i]

            bg_estimates.append(bg)
        
        self.peak_finder.newImages(images)
        self.peak_fitter.newImages(images)

        peaks = False
        for i in range(self.peak_finder.iterations):

            # Update background estimate.
            for j in range(self.n_planes):
                self.peak_finder.subtractBackground(images[j] - fit_peaks_images[j], bg_estimates[j], j)

            # Find new peaks.
            [found_new_peaks, peaks] = self.peak_finder.findPeaks(fit_peaks_images, peaks)

            # Fit new peaks.
            if isinstance(peaks, numpy.ndarray):
                [peaks, fit_peaks_images] = self.peak_fitter.fitPeaks(peaks)

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
        #peaks[:,utilC.getStatusIndex()] = 1.0
        converged_peaks = super().getConvergedPeaks(peaks)
        return self.peak_fitter.rescaleZ(converged_peaks)
