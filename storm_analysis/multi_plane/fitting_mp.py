#!/usr/bin/env python
"""
(Standard) multi-plane peak finder and fitter. This is used for both arbitrary PSF
and 3D-DAOSTORM Gaussian PSF fitting.

This provides sub-classes of the classes in storm_analysis.sa_library.fitting

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

import storm_analysis.multi_plane.mp_utilities as mpUtil

import storm_analysis.sa_library.affine_transform_c as affineTransformC
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC


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
    def __init__(self, parameters = None, **kwds):
        kwds["parameters"] = parameters
        super(MPPeakFinder, self).__init__(**kwds)

        self.atrans = [None]
        self.backgrounds = []
        self.bg_filters = []           # Matched filters for background estimation.
        self.images = []
        self.mfilters = []             # Matched filters for foreground estimation.
        self.n_channels = None
        self.variances = []            # sCMOS variance data.
        self.vfilters = []             # Matched filters for (background) variance estimation.
        self.xt = []
        self.yt = []
        self.z_values = []
    
    def backgroundEstimator(self, image, channel):
        return self.bg_filters[channel].convolve(image)

    def cleanUp(self):
        super(MPPeakFinder, self).cleanUp()

        # Clean up transforms.
        for at in self.atrans:
            if at is not None:
                at.cleanup()

        # Clean up foreground filters and variance filters.
        for i in range(len(self.mfilters)):
            for j in range(len(self.mfilters[i])):
                self.mfilters[i][j].cleanup()
                self.vfilters[i][j].cleanup()

        # Clean up background filters.
        for i in range(len(self.bg_filters)):
            self.bg_filters[i].cleanup()

    def estimateBackground(self, fit_peaks_images, bg_estimates):
        """
        Estimate the background for the images.
        """
        if bg_estimates[0] is not None:
            self.backgrounds = bg_estimates

        else:
            self.backgrounds = []
            for i in range(len(self.images)):
                self.backgrounds.append(self.backgroundEstimator(self.images[i] - fit_peaks_images[i], i))

        # Save results if requested.
        if self.check_mode:
            with tifffile.TiffWriter("bg_estimate.tif") as tf:
                for bg in self.backgrounds:
                    tf.save(bg.astype(numpy.float32))

        return self.backgrounds
            
    def loadMapping(self):
        """
        Load the channel to channel mapping and create affine transform objects.
        """
        mappings = {}
        if self.parameters.hasAttr("mapping"):
            if os.path.exists(self.parameters.getAttr("mapping")):
                with open(self.parameters.getAttr("mapping"), 'rb') as fp:
                    mappings = pickle.load(fp)
            else:
                raise Exception("Mapping file " + self.parameters.getAttr("mapping") + " does not exist.")

        for i in range(self.n_channels-1):
            self.xt.append(mpUtil.marginCorrect(mappings["0_" + str(i+1) + "_x"], self.margin))
            self.yt.append(mpUtil.marginCorrect(mappings["0_" + str(i+1) + "_y"], self.margin))
            self.atrans.append(affineTransformC.AffineTransform(xt = self.xt[i],
                                                                yt = self.yt[i]))
            
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
                    tf.save(at_image.astype(numpy.float32))

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
                for fi in fit_peaks_images:
                    tf.save(fi.astype(numpy.float32))

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

        # Remove problematic values.
        for bg in bg_variances:
            mask = (bg <= 0.1)
            if (numpy.sum(mask) > 0):
                if self.check_mode:
                    print("Warning! small and/or negative values detected in background variance.")
                bg[mask] = 0.1
        
        # Save results if needed for debugging purposes.
        if self.check_mode:
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
                    tf.save(fg.astype(numpy.float32))

            with tifffile.TiffWriter("fg_bg_ratio.tif") as tf:
                for fg_bg_ratio in fg_bg_ratios:
                    tf.save(fg_bg_ratio.astype(numpy.float32))

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
        """
        We initialize the following here because at __init__ we
        don't know how big the images are. This is called after
        the analysis specific version pads the variances and
        creates the mfilters[] and the vfilters[] class members.
        
        Note the assumption that every frame in all the movies
        is the same size.
        """
        
        # Create mask to limit peak finding to a user defined sub-region of the image.
        #
        self.peak_mask = fitting.peakMask(variances[0].shape, self.parameters, self.margin)

        # Create matched filter for background. There is one of these for
        # each imaging plane for the benefit of memoization.
        #
        for i in range(self.n_channels):
            bg_psf = fitting.gaussianPSF(variances[0].shape, self.parameters.getAttr("background_sigma"))
            self.bg_filters.append(matchedFilterC.MatchedFilter(bg_psf,
                                                                fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                                memoize = True,
                                                                max_diff = 1.0e-3))

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
        #
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
                    tf.save(var.astype(numpy.float32))

 
    

class MPPeakFinderArb(MPPeakFinder):
    """
    Multi-plane peak finding for the spline, pupil function or PSF FFT fitting models.
    """
    def __init__(self, parameters = None, psf_objects = None, **kwds):
        kwds["parameters"] = parameters
        super(MPPeakFinderArb, self).__init__(**kwds)

        self.mfilters_z = []
        self.n_channels = len(psf_objects)
        self.psf_objects = psf_objects

        # Assert that all the PSF objects are the same size.
        #
        for i in range(1, len(self.psf_objects)):
            assert(self.psf_objects[0].getSize() == self.psf_objects[i].getSize())

        # Update margin based on the psf object size.
        #
        self.margin = self.psf_objects[0].getMargin()

        # Note the assumption that the PSFs for each plane all use
        # the same z scale / have the same z range.
        #
        # 'z_value' is in units of microns.
        # self.mfilters_z is the Z position in nanometers.
        #
        self.mfilters_z = list(map(lambda x: x * 1.0e+3, parameters.getAttr("z_value", [0.0])))
        for zval in self.mfilters_z:
            assert self.psf_objects[0].isValidZ(zval)
            self.z_values.append(self.psf_objects[0].getScaledZ(zval))

        # Configure maxima finder.
        #
        self.mfinder = iaUtilsC.MaximaFinder(margin = self.margin,
                                             radius = self.find_max_radius,
                                             threshold = self.threshold,
                                             z_values = self.z_values)
            
        # Load the channel to channel mapping file. We need the correct value
        # for self.margin for this.
        #
        self.loadMapping()

        # Load pre-specified peak locations, if any.
        #
        if parameters.hasAttr("peak_locations"):
            [self.peak_locations, self.peak_locations_type] = fitting.getPeakLocations(parameters.getAttr("peak_locations"),
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
                self.peak_locations["z"] = self.psf_objects[0].getScaledZ(self.peak_locations["z"])

        # A lot additional initialization occurs in the setVariances method
        # as this is the first time we'll know the image size.
        #

    def setVariances(self, variances):
        """
        setVariances() customized for arbitrary PSFs.
        """

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

        # Create "foreground" and "variance" filters.
        #
        # These are stored in a list indexed by z value, then by
        # channel / plane. So self.mfilters[1][2] is the filter
        # for z value 1, plane 2.
        #
        for i, mfilter_z in enumerate(self.mfilters_z):
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
                self.mfilters[i].append(matchedFilterC.MatchedFilter(psf_norm,
                                                                     fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                                     memoize = True,
                                                                     max_diff = 1.0e-3))
                self.vfilters[i].append(matchedFilterC.MatchedFilter(psf_norm * psf_norm,
                                                                     fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                                     memoize = True,
                                                                     max_diff = 1.0e-3))

                # Save a pictures of the PSFs for debugging purposes.
                if self.check_mode:
                    print("psf max", numpy.max(psf))
                    filename = "psf_z{0:.3f}_c{1:d}.tif".format(mfilter_z, j)
                    tifffile.imsave(filename, psf.astype(numpy.float32))

        # This handles the rest of the initialization.
        #
        super(MPPeakFinderArb, self).setVariances(variances)

        return variances


class MPPeakFinderDao(MPPeakFinder):
    """
    Multi-plane peak finding for 3D-DAOSTORM (Gaussian) PSF fitting model.

    In theory this could follow MPPeakFinderArb() and do convolutions at 
    different Z values to determine a better guess for the peak starting
    sigmas. For now though this is being kept simple and it just does 
    the convolution at a single plane and uses the 'foreground_sigma'
    parameter for the gaussian sigma.
    """
    def __init__(self, parameters = None, n_channels = None, **kwds):
        kwds["parameters"] = parameters
        super(MPPeakFinderDao, self).__init__(**kwds)

        self.n_channels = n_channels
        self.z_values = [0.0]
        
        # Figure out what margin and ROI to use.
        if (self.parameters.getAttr("roi_size", -1) != -1):
            self.roi_size = parameters.getAttr("roi_size")
        else:
            
            # Calculate roi size based on sigma.
            self.roi_size = int(20.0 * self.sigma)

        self.margin = int(self.roi_size/2 + 2)

        # Configure maxima finder.
        #
        self.mfinder = iaUtilsC.MaximaFinder(margin = self.margin,
                                             radius = self.find_max_radius,
                                             threshold = self.threshold,
                                             z_values = self.z_values)
            
        # Load the channel to channel mapping file. We need the correct value
        # for self.margin for this.
        #
        self.loadMapping()

        # Load pre-specified peak locations, if any.
        #
        if parameters.hasAttr("peak_locations"):
            [self.peak_locations, self.peak_locations_type] = fitting.getPeakLocations(parameters.getAttr("peak_locations"),
                                                                                       self.margin,
                                                                                       parameters.getAttr("pixel_size"),
                                                                                       self.sigma)

            # Set initial z value to 0.0. This is just a placeholder, this
            # value won't be used for anything.
            #
            self.peak_locations["z"][:] = self.z_values[0]

        # A lot additional initialization occurs in the setVariances method
        # as this is the first time we'll know the image size.
        #

    def peakFinder(self, fit_peaks_image):
        # Add sigma field to the peak initial parameters.
        #
        peaks = super(MPPeakFinderDao, self).peakFinder(fit_peaks_image)
        peaks["sigma"] = numpy.ones(peaks["x"].size) * self.sigma
        return peaks

    def setVariances(self, variances):
        """
        setVariances() customized for gaussian PSFs.
        """

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

        # Create "foreground" and "variance" filters. There is
        # only one z value here.
        #
        # These are stored in a list indexed by z value, then by
        # channel / plane. So self.mfilters[1][2] is the filter
        # for z value 1, plane 2.
        #
        self.mfilters.append([])
        self.vfilters.append([])

        psf_norm = fitting.gaussianPSF(variances[0].shape, self.parameters.getAttr("foreground_sigma"))
        var_norm = psf_norm * psf_norm

        for i in range(self.n_channels):
            self.mfilters[0].append(matchedFilterC.MatchedFilter(psf_norm,
                                                                 fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                                 memoize = True,
                                                                 max_diff = 1.0e-3))
            self.vfilters[0].append(matchedFilterC.MatchedFilter(var_norm,
                                                                 fftw_estimate = self.parameters.getAttr("fftw_estimate"),
                                                                 memoize = True,
                                                                 max_diff = 1.0e-3))

            # Save a pictures of the PSFs for debugging purposes.
            if self.check_mode:
                print("psf max", numpy.max(psf))
                filename = "psf_z0.0_c{1:d}.tif".format(j)
                tifffile.imsave(filename, psf.astype(numpy.float32))

        # This handles the rest of the initialization.
        #
        super(MPPeakFinderDao, self).setVariances(variances)

        return variances

    
class MPPeakFitterArb(fitting.PeakFitterArbitraryPSF):
    """
    Multi-plane peak fitting, arbitrary PSF.
    """
    def getPeakProperty(self, pname, channel = 0):
        return self.mfitter.getPeakProperty(pname, channel = channel)


class MPPeakFitterDao(fitting.PeakFitter):
    """
    Multi-plane peak fitting, guassian PSF.
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
            temp = {}
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

    
