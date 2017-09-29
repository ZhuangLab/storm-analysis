#!/usr/bin/env python
"""
Cubic spline peak finder.

Hazen 03/16
"""

import pickle
import numpy

import tifffile

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC

import storm_analysis.spliner.cubic_fit_c as cubicFitC
import storm_analysis.spliner.spline_to_psf as splineToPSF


class SplinerPeakFinder(fitting.PeakFinder):
    """
    Spliner peak finding.
    """
    def __init__(self, parameters = None, **kwds):
        kwds["parameters"] = parameters
        super(SplinerPeakFinder, self).__init__(self, parameters)

        self.height_rescale = []
        self.fg_mfilter = []
        self.fg_mfilter_zval = []
        self.fg_vfilter = []
        self.z_values = []

        # Load the spline.
        self.s_to_psf = splineToPSF.loadSpline(parameters.getAttr("spline"))

        # Update margin based on the spline size.
        old_margin = self.margin
        self.margin = int((self.s_to_psf.getSize() + 1)/4 + 2)

        self.fg_mfilter_zval = parameters.getAttr("z_value", [0.0])
        for zval in self.fg_mfilter_zval:
            self.z_values.append(self.s_to_psf.getScaledZ(zval))

        if parameters.hasAttr("peak_locations"):

            # Correct for any difference in the margins.
            self.peak_locations[:,utilC.getXCenterIndex()] += self.margin - old_margin
            self.peak_locations[:,utilC.getYCenterIndex()] += self.margin - old_margin

            # Provide the "correct" starting z value.
            self.peak_locations[:,utilC.getZCenterIndex()] = self.z_value[0]

    def newImage(self, new_image):
        """
        This is called once at the start of the analysis of a new image.
        
        new_image - A 2D numpy array.
        """
        super(SplinerPeakFinder, self).newImage(new_image)
        
        #
        # If does not already exist, create filter objects from
        # the best fit spline at different z value of the PSF.
        #
        # As not all PSFs will be maximal in the center we can't just
        # use the image intensity at the center as the starting
        # height value. Instead we will use the intensity at the
        # peak center of the convolved image, then adjust this
        # value by the height_rescale parameter.
        #
        if (len(self.fg_mfilter) == 0):
            for zval in self.fg_mfilter_zval:
                psf = self.s_to_psf.getPSF(zval,
                                           shape = new_image.shape,
                                           normalize = False)
                psf_norm = psf/numpy.sum(psf)                
                self.fg_mfilter.append(matchedFilterC.MatchedFilter(psf_norm))
                self.fg_vfilter.append(matchedFilterC.MatchedFilter(psf_norm * psf_norm))

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
                self.height_rescale.append(1.0/numpy.sum(psf * psf_norm))

                # Save a picture of the PSF for debugging purposes.
                if self.check_mode:
                    print("psf max", numpy.max(psf))
                    temp = 10000.0 * psf + 100.0
                    filename = "psf_{0:.3f}.tif".format(zval)
                    tifffile.imsave(filename, temp.astype(numpy.uint16))

        self.taken = []
        for i in range(len(self.mfilter)):
            self.taken.append(numpy.zeros(new_image.shape, dtype=numpy.int32))
            
    def peakFinder(self, fit_peaks_image):
        """
        This method does the actual peak finding.
        """
        all_new_peaks = None

        # Calculate background variance.
        #
        # Note the assumption here that we are working in units of photo-electrons
        # so Poisson statistics applies, variance = mean.
        #
        bg_var = self.background + fit_peaks_image

        # Add camera variance if set.
        if self.camera_variance is not None:
            bg_var += self.camera_variance
            
        #
        # Find peaks in image convolved with the PSF at different z values.
        #
        for i in range(len(self.fg_mfilter)):

            # Estimate background variance at this particular z value.
            bg_var = self.fg_vfilter[i].convolve(bg_var)

            # Check for problematic values.
            #
            # Note: numpy will also complain when we try to take the sqrt of a negative number.
            #
            if self.check_mode:            
                mask = (bg_var <= 0.0)
                if (numpy.sum(mask) > 0):
                    print("Warning! zero and/or negative values detected in background variance!")
                    
            # Convert to standard deviation.
            bg_std = numpy.sqrt(bg_var)

            # Calculate foreground.
            foreground = self.image - self.background - fit_peaks_image
            foreground = self.fg_mfilter[i].convolve(foreground)

            if self.check_mode:
                with tifffile.TiffWriter("foreground_{0:.2f}.tif".format(self.z_values[i])) as tf:
                    tf.save(numpy.transpose(foreground.astype(numpy.float32)))
                    
            # Calculate foreground in units of signal to noise.
            fg_bg_ratio = foreground/bg_std
        
            if self.check_mode:
                with tifffile.TiffWriter("fg_bg_ratio_{0:.2f}.tif".format(self.z_values[i])) as tf:
                    tf.save(numpy.transpose(fg_bg_ratio.astype(numpy.float32)))        
                    
            # Mask the image so that peaks are only found in the AOI.
            masked_image = fg_bg_ratio * self.peak_mask
        
            # Identify local maxima in the masked ratio image.
            [new_peaks, taken] = utilC.findLocalMaxima(fg_bg_ratio,
                                                       taken,
                                                       self.cur_threshold,
                                                       self.find_max_radius,
                                                       self.margin)

            # Fill in initial values for peak height, background and sigma.
            new_peaks = utilC.initializePeaks(new_peaks,                    # The new peaks.
                                              foreground + self.background, # Convolved image + background.
                                              self.background,              # The current estimate of the background.
                                              self.sigma,                   # The starting sigma value.
                                              z_value)                      # The starting z value.

            # Correct initial peak heights, self.height_rescale is an estimate
            # of the effect of PSF convolution on the height of the original
            # localization.
            h_index = utilC.getHeightIndex()
            new_peaks[:,h_index] = new_peaks[:,h_index] * self.height_rescale[i]
            
            if all_new_peaks is None:
                all_new_peaks = new_peaks
            else:
                all_new_peaks = numpy.append(all_new_peaks, new_peaks, axis = 0)
                
        #
        # Remove the dimmer of two peaks with similar x,y values but different z values.
        #
        if (len(self.mfilter) > 1):

            if self.check_mode:
                print("Before peak removal", all_new_peaks.shape)
                for i in range(all_new_peaks.shape[0]):
                    print(all_new_peaks[i,:])
                print("")
            
            all_new_peaks = utilC.removeClosePeaks(all_new_peaks,                                               
                                                   self.find_max_radius,
                                                   self.find_max_radius)

            if self.check_mode:
                print("After peak removal", all_new_peaks.shape)
                for i in range(all_new_peaks.shape[0]):
                    print(all_new_peaks[i,:])
                print("")

        return all_new_peaks


class SplinerPeakFitter(fitting.PeakFitter):
    """
    Spliner peak fitting.
    """
    # Convert from spline z units to real z units.
    def rescaleZ(self, peaks):
        return self.mfitter.rescaleZ(peaks)
        

class SplinerFinderFitter(fitting.PeakFinderFitter):
    """
    Class for spline based peak finding and fitting.
    """
    def __init__(self, **kwds):
        super(SplinerFinderFitter, self).__init__(**kwds)

        # Update margin.
        self.margin = self.peak_finder.margin

    def getConvergedPeaks(self, peaks):
        converged_peaks = fitting.PeakFinderFitter.getConvergedPeaks(self, peaks)
        return self.peak_fitter.rescaleZ(converged_peaks)


def initFindAndFit(parameters):
    """
    Initialize and return a SplinerFinderFitter object.
    """
    # Create peak finder.
    finder = SplinerPeakFinder(parameters = parameters)

    # Load variance, scale by gain.
    #
    # Variance is in units of ADU*ADU.
    # Gain is ADU/photo-electron.
    #
    variance = None
    if parameter.hasAttr("camera_calibration"):
        [offset, variance, gain] = numpy.load(parameters.getAttr("camera_calibration"))
        variance = variance/(gain*gain)

        # Set variance in the peak finder, this method also pads the
        # variance to the correct size.
        variance = finder.setVariance(variance)

    # Load spline and create the appropriate type of spline fitter.
    with open(parameters.getAttr("spline"), 'rb') as fp:
        psf_data = pickle.load(fp)

    # If this is a 3D spline, get the range it covers in microns.
    if(psf_data["type"] == "3D"):
        min_z = psf_data["zmin"]/1000.0
        max_z = psf_data["zmax"]/1000.0
            
    spline = psf_data["spline"]
    coeff = psf_data["coeff"]

    # Create C fitter object.
    if (len(spline.shape) == 2):
        mfitter = cubicFitC.CSpline2DFit(scmos_cal = variance,
                                         spline_vals = spline,
                                         coeff_vals = coeff)
    else:
        mfitter = cubicFitC.CSpline3DFit(scmos_cal = variance,
                                         spline_vals = spline,
                                         coeff_vals = coeff,
                                         min_z = min_z,
                                         max_z = max_z)

    # Create peak fitter.
    fitter = SplinerPeakFitter(mfitter = mfitter,
                               parameters = parameters)

    return SplinerFinderFitter(peak_finder = finder,
                               peak_fitter = fitter)    
