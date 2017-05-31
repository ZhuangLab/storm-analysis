#!/usr/bin/env python
"""
(Standard) multi-plane peak finder and fitter.

"channel", "frame", "image" and "plane" are used somewhat interchangeably..

Hazen 05/17
"""

import pickle
import numpy
import tifffile

import storm_analysis.sa_library.affine_transform_c as affineTransformC
import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC

import storm_analysis.sa_utilities.std_analysis as stdAnalysis

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
        return self.backgrounds[plane]

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

    All locations are relative to plane 0.

    The expectation is that we are working with an image that has already 
    been corrected for gain and offset.

    This works with affine transformed versions of the images so that they
    can all be overlaid on top of each other.
    """
    def __init__(self, parameters):
        super().__init__(parameters)

        self.atrans = [None]
        self.backgrounds = []
        self.height_rescale = []
        self.images = []
        self.mfilters = []
        self.mfilters_z = []
        self.n_channels = 0
        self.s_to_psfs = []
        self.variances = []
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
        with open(parameters.getAttr("mapping"), 'rb') as fp:
            mappings = pickle.load(fp)

        for i in range(self.n_channels-1):
            xt = mappings["0_" + str(i+1) + "_x"]
            yt = mappings["0_" + str(i+1) + "_y"]
            self.atrans.append(affineTransformC.AffineTransform(margin = self.margin,
                                                                xt = xt,
                                                                yt = yt))

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

    def cleanUp(self):
        for at in self.atrans:
            if at is not None:
                at.cleanup()
        
    def findPeak(self, no_bg_image, peaks):
        pass
    
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
        for i in range(len(self.mfilters)):
            self.taken.append(numpy.zeros(new_images[0].shape, dtype=numpy.int32))

        #
        # Apply affine transforms to input images.
        #
        self.images = []
        for i in range(self.n_channels):
            if self.atrans[i] is None:
                self.images.append(new_images[i].copy())
            else:
                self.images.append(self.atrans[i].transform(new_images[i]))

        # For checking that we're doing the transform correctly.
        if True:
            with tifffile.TiffWriter("transform.tif") as tf:
                for image in self.images:
                    tf.save(image.astype(numpy.float32))
        
        #
        # We initialize the following here because at __init__ we
        # don't know how big the images are.
        #
        
        # Create mask to limit peak finding to a user defined sub-region of the image.
        if self.peak_mask is None:
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
        # This is basically the same as spliner, except that we
        # need filters for each z value and for each plane.
        #
        if (len(self.mfilters) == 0):
            for i, mfilter_z in enumerate(self.mfilters_z):
                self.mfilters.append([])
                h_rescale = 0.0
                for j, s_to_psf in enumerate(self.s_to_psfs):
                    psf = s_to_psf.getPSF(mfilter_z,
                                          shape = new_images[0].shape,
                                          normalize = False)
                    psf_norm = psf/numpy.sum(psf)
                    h_rescale += numpy.sum(psf * psf_norm)
                    self.mfilters[i].append(matchedFilterC.MatchedFilter(psf_norm))

                    # Save a pictures of the PSFs for debugging purposes.
                    if False:
                        print("psf max", numpy.max(psf))
                        temp = 10000.0 * psf + 100.0
                        filename = "psf_z{0:.3f}_c{1:d}.tif".format(mfilter_z, j)
                        tifffile.imsave(filename, temp.astype(numpy.uint16))
                    
                self.height_rescale.append(1.0/h_rescale)
                
    def peakFinder(self, no_bg_image):
        pass

    def setVariances(self, variances):
        self.variances = variances

    def subtractBackground(self, image, bg_estimate):
        pass


class MPPeakFitter(fitting.PeakFitter):
    """
    Multi-plane peak fitting.
    """
    def __init__(self, parameters):
        super().__init__(parameters)
        self.mp_fitter = None
        self.variances = []
        
        # Load the plane to plane mapping.
        
        # Load the splines & create the multi-plane fitter.
        self.coeffs = []
        self.splines = []
        for spline_name in getSplineAttrs(parameters):
            with open(parameters.getAttr(spline_name), 'rb') as fp:
                spline_data = pickle.load(fp)
            self.zmin = spline_data["zmin"]/1000.0
            self.zmax = spline_data["zmax"]/1000.0
            self.coeffs.append(spline_data["coeff"])
            self.splines.append(spline_data["spline"])

    def cleanUp(self):
        pass
    
    def fitPeaks(self, peaks):
        pass

    def newImages(self, new_images):
        pass

    def rescaleZ(self, peaks):
        """
        Convert from spline z units to real z units.
        """
        return self.mp_fitter.rescaleZ(peaks, self.zmin, self.zmax)

    def setVariances(self, variances):
        self.variances = variances


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

    def analyzeImage(self, movie_reader, save_residual = False, verbose = False):

        # Load and scale the images.
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

        residual = numpy.zeros(images[0].shape)
        
        self.peak_finder.newImages(images)
        self.peak_fitter.newImages(images)

        peaks = False
        for i in range(self.peak_finder.iterations):
            continue

        if isinstance(peaks, numpy.ndarray):
            peaks[:,utilC.getXCenterIndex()] -= float(self.margin)
            peaks[:,utilC.getYCenterIndex()] -= float(self.margin)

        return [peaks, residual]

    def cleanUp(self):
        self.peak_finder.cleanUp()
        self.peak_fitter.cleanUp()

    def getConvergedPeaks(self, peaks):
        pass
#        converged_peaks = super().getConvergedPeaks(self, peaks)
#        return self.peak_fitter.rescaleZ(converged_peaks)
