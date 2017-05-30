#!/usr/bin/env python
"""
Multi-plane peak finder and fitter.

Hazen 05/17
"""

import pickle
import numpy

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
        self.planes = []
        
        # Load the movies and offsets for each plane/channel.
        for ext in getExtAttrs(parameters):
            self.planes.append(datareader.inferReader(base_name + ext))

        for offset in getOffsetAttrs(parameters):
            self.offsets.append(parameters.getAttr("channel" + str(i) + "_offset"))

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
            for i, bg_estimator in enumerate(self.bg_estimators):
                self.backgrounds[i] = bg_estimator.estimateBG(self.cur_frame + self.offsets[i])

            # Load planes. Zero / negative value removal is done later in the
            # analysis (after offset subtraction and gain correction).
            for i, plane in enumerate(self.planes):
                self.frames[i] = plane.loadAFrame(self.cur_frame + self.offsets[i])

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
                    self.cur_frame = parameters.getAttr("start_frame")
        
        # Figure out where to stop.
        self.max_frame = self.movie_l
        if self.parameters.hasAttr("max_frame"):
            if (self.parameters.getAttr("max_frame") > 0):
                if (self.parameters.getAttr("max_frame") < self.movie_l):
                    self.max_frame = parameters.getAttr("max_frame")

        # Configure background estimator, if any.
        if (self.parameters.getAttr("static_background_estimate", 0) > 0):
            print("Using static background estimator.")
            s_size = self.parameters.getAttr("static_background_estimate")
            for i in range(len(self.planes)):
                bg_est = static_background.StaticBGEstimator(self.planes[i],
                                                             start_frame = self.cur_frame + self.offsets[i],
                                                             sample_size = s_size))
                self.bg_estimators.append(bg_est)


class MPPeakFinder(fitting.PeakFinder):
    """
    Multi-plane peak finding.

    All locations are relative to plane 0.

    Also, the expectation is that we are working with an image 
    that has already been corrected for gain and offset.
    """
    def __init__(self, parameters):
        super().__init__(parameters)

        self.height_rescale = []
        self.mfilters = []
        self.mfilters_z = []
        self.s_to_psfs = []
        self.z_values = []

        # Load the plane to plane mapping data.

        # Load the splines.
        for spline_name in getSplineAttrs(parameters):
            self.s_to_psfs.append(splineToPSF.loadSpine(parameters.getAttr(spline_name)))

        # Update margin based on the spline size. Note the assumption
        # that all of the splines are the same size, or at least smaller
        # than spline for plane 0.
        old_margin = self.margin
        self.margin = int((self.s_to_psfs[0].getSize() + 1)/4 + 2)

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

    def findPeak(self, no_bg_image, peaks):
        pass
    
    def newImage(self, new_image):
        super().newImage(new_image)

        #
        # This is basically the same as spliner, except that we
        # need filters for each z value and for each plane.
        #
        if (len(self.mfilters) == 0):
            for i, mfilter_z in enumerate(self.mfilter_z):
                self.mfilters.append([])
                h_rescale = 0.0
                for j, s_to_psf in enumerate(self.s_to_psfs):
                    psf = s_to_psf.getPSF(mfilter_z,
                                          shape = new_image.shape,
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

        self.taken = []
        for i in range(len(self.mfilters)):
            self.taken.append(numpy.zeros(new_image.shape, dtype=numpy.int32))

    def peakFinder(self, no_bg_image):
        pass

    def subtractBackground(self, image, bg_estimate):
        pass


class MPPeakFitter(fitting.PeakFitter):
    """
    Multi-plane peak fitting.
    """
    def __init__(self, parameters):
        super().__init__(parameters)
        self.mp_fitter = None
        
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

    def fitPeaks(self, peaks):
        pass

    def newImage(self, image):
        pass

    def rescaleZ(self, peaks):
        """
        Convert from spline z units to real z units.
        """
        return self.mp_fitter.rescaleZ(peaks, self.zmin, self.zmax)


class MPFinderFitter(fitting.PeakFinderFitter):
    """
    Multi-plane peak finding and fitting.
    """
    def __init__(self, parameters):
        super().__init__(parameters)
        self.peak_finder = MPPeakFinder(parameters)
        self.peak_fitter = MPPeakFitter(parameters)

        # Update margin.
        self.margin = self.peak_finder.margin

        # Load sCMOS calibration data.

    def analyzeImage(self, movie_reader, save_residual = False, verbose = False):
        pass

    def getConvergedPeaks(self, peaks):
        pass
#        converged_peaks = super().getConvergedPeaks(self, peaks)
#        return self.peak_fitter.rescaleZ(converged_peaks)
