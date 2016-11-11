#!/usr/bin/python
#
# Cubic spline peak finder.
#
# Hazen 03/16
#

import pickle
import numpy

import tifffile

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC

import storm_analysis.spliner.cubic_fit_c as cubicFitC
import storm_analysis.spliner.spline_to_psf as splineToPSF


#
# Spliner peak finding.
#
class SplinerPeakFinder(fitting.PeakFinder):

    def __init__(self, parameters):
        fitting.PeakFinder.__init__(self, parameters)
        self.height_rescale = []
        self.mfilter = []
        self.mfilter_z = []
        self.z_value = []
        self.s_to_psf = splineToPSF.SplineToPSF(parameters.spline)

        # Update margin based on the spline size.
        old_margin = self.margin
        self.margin = int((self.s_to_psf.getSize() + 1)/4 + 2)
        
        if hasattr(parameters, "z_value"):

            #
            # If your PSF is not degenerate* in Z then it could be helpful
            # to try multiple z starting values. However most common 3D
            # PSFs such as astigmatism do not meet this criteria. The only
            # one that does meet this criteria that is in (sort of) common
            # use is the double-helix PSF.
            #
            # * By degenerate I mean that the PSF at one z value can
            #   be modeled (with reasonable accuracy) by summing several
            #   PSFs with a different z value. For example, most
            #   astigmatic PSFs z != 0 can be modeled by
            #   summing several z = 0 PSFs with variable x,y positions.
            #
            if isinstance(parameters.z_value, list):
                for z_value in parameters.z_value:
                    self.mfilter_z.append(z_value)
                    self.z_value.append(self.s_to_psf.getScaledZ(z_value))
            else:
                self.mfilter_z = [parameters.z_value]
                self.z_value = [self.s_to_psf.getScaledZ(parameters.z_value)]
        else:
            self.mfilter_z = [0.0]
            self.z_value = [self.s_to_psf.getScaledZ(0.0)]

        if hasattr(parameters, "peak_locations"):

            # Correct for any difference in the margins.
            self.peak_locations[:,utilC.getXCenterIndex()] += self.margin - old_margin
            self.peak_locations[:,utilC.getYCenterIndex()] += self.margin - old_margin

            # Provide the "correct" starting z value.
            self.peak_locations[:,utilC.getZCenterIndex()] = self.z_value[0]

    def newImage(self, new_image):
        fitting.PeakFinder.newImage(self, new_image)
        
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
        if (len(self.mfilter) == 0):
            for mfilter_z in self.mfilter_z:
                psf = self.s_to_psf.getPSF(mfilter_z,
                                           shape = new_image.shape,
                                           normalize = False)
                psf_norm = psf/numpy.sum(psf)
                self.height_rescale.append(1.0/numpy.sum(psf * psf_norm))
                self.mfilter.append(matchedFilterC.MatchedFilter(psf_norm))

                # Save a picture of the PSF for debugging purposes.
                if False:
                    print("psf max", numpy.max(psf))
                    temp = 10000.0 * psf + 100.0
                    filename = "psf_{0:.3f}.tif".format(mfilter_z)
                    tifffile.imsave(filename, temp.astype(numpy.uint16))

        self.taken = []
        for i in range(len(self.mfilter)):
            self.taken.append(numpy.zeros(new_image.shape, dtype=numpy.int32))
                           
    def peakFinder(self, no_bg_image):

        all_new_peaks = None

        save_convolution = True
        print("making image")
        if save_convolution:
            tif = tifffile.TiffWriter("image.tif")
            tif.save(self.image.astype(numpy.uint16) + 100)
            tif.save(no_bg_image.astype(numpy.uint16) + 100)
            tif.save(self.background.astype(numpy.uint16) + 100)
            
        #
        # Find peaks in image convolved with the PSF at different z values.
        #
        for i in range(len(self.mfilter)):

            height_rescale = self.height_rescale[i]
            mfilter = self.mfilter[i]
            taken = self.taken[i]
            z_value = self.z_value[i]

            # Smooth image with gaussian filter.
            smooth_image = mfilter.convolve(no_bg_image)

            if save_convolution:
                tif.save(smooth_image.astype(numpy.uint16) + 100)

            # Mask the image so that peaks are only found in the AOI.
            masked_image = smooth_image * self.peak_mask
        
            # Identify local maxima in the masked image.
            [new_peaks, taken] = utilC.findLocalMaxima(masked_image,
                                                       taken,
                                                       self.cur_threshold,
                                                       self.find_max_radius,
                                                       self.margin)

            #
            # Fill in initial values for peak height, background and sigma.
            #
            # FIXME: We just add the smoothed image and the background together
            #        as a hack so that we can still use the initializePeaks()
            #        function.
            #
            new_peaks = utilC.initializePeaks(new_peaks,                      # The new peaks.
                                              smooth_image + self.background, # Smooth image + background.
                                              self.background,                # The current estimate of the background.
                                              self.sigma,                     # The starting sigma value.
                                              z_value)                        # The starting z value.

            # Correct initial peak heights.
            h_index = utilC.getHeightIndex()
            new_peaks[:,h_index] = new_peaks[:,h_index] * height_rescale
            
            if all_new_peaks is None:
                all_new_peaks = new_peaks
            else:
                all_new_peaks = numpy.append(all_new_peaks, new_peaks, axis = 0)

        if save_convolution:
            tif.close()
                
        #
        # Remove the dimmer of two peaks with similar x,y values but different z values.
        #
        if (len(self.mfilter) > 1):

            if False:
                print("before", all_new_peaks.shape)
                for i in range(all_new_peaks.shape[0]):
                    print(all_new_peaks[i,:])
                print("")
            
            all_new_peaks = utilC.removeClosePeaks(all_new_peaks,                                               
                                                   self.find_max_radius,
                                                   self.find_max_radius)

            if False:
                print("after", all_new_peaks.shape)
                for i in range(all_new_peaks.shape[0]):
                    print(all_new_peaks[i,:])
                print("")
                
        return all_new_peaks


#
# Spliner peak fitting.
#
class SplinerPeakFitter(fitting.PeakFitter):

    def __init__(self, parameters):
        fitting.PeakFitter.__init__(self, parameters)

        # Load spline and create the appropriate type of spline fitter.
        psf_data = pickle.load(open(parameters.spline, 'rb'))
        self.zmin = psf_data["zmin"]/1000.0
        self.zmax = psf_data["zmax"]/1000.0
        self.spline = psf_data["spline"]

        save_coeff = True
        self.coeff = False
        if ("coeff" in psf_data):
            save_coeff = False
            self.coeff = psf_data["coeff"]

        if (len(self.spline.shape)==2):
            self.spline_type = "2D"
            self.mfitter = cubicFitC.CSpline2DFit(self.spline, self.coeff, None)
        else:
            self.spline_type = "3D"
            self.mfitter = cubicFitC.CSpline3DFit(self.spline, self.coeff, None)

        # Save the coefficients for faster start up.
        if save_coeff:
            psf_data["coeff"] = self.sfitter.getCoeff()
            pickle.dump(psf_data, open(parameters.spline, "w"))

    def fitPeaks(self, peaks):
        
        # Fit to update peak locations.
        [fit_peaks, residual] = self.peakFitter(peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks,
                                              min_height = 0.9*self.threshold)

        # Remove peaks that are too close to each other & refit.
        fit_peaks = utilC.removeClosePeaks(fit_peaks, self.sigma, self.neighborhood)
        [fit_peaks, residual] = self.peakFitter(fit_peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks,
                                              min_height = 0.9*self.threshold)

        return [fit_peaks, residual]

    def newImage(self, image):
        self.mfitter.newImage(image)

    # Convert from spline z units to real z units.
    def rescaleZ(self, peaks):
        if (self.spline_type == "3D"):
            return self.mfitter.rescaleZ(peaks, self.zmin, self.zmax)
        else:
            return peaks

#
# Class to encapsulate spline based peak finding and fitting.
#
class SplinerFinderFitter(fitting.PeakFinderFitter):

    def __init__(self, parameters):
        fitting.PeakFinderFitter.__init__(self, parameters)
        self.peak_finder = SplinerPeakFinder(parameters)
        self.peak_fitter = SplinerPeakFitter(parameters)

        # Update margin.
        self.margin = self.peak_finder.margin

    def analyzeImage(self, new_image, bg_estimate = None, save_residual = False, verbose = False):
        return fitting.PeakFinderFitter.analyzeImage(self,
                                                     new_image,
                                                     bg_estimate = bg_estimate,
                                                     save_residual = save_residual,
                                                     verbose = verbose)

    def getConvergedPeaks(self, peaks):
        converged_peaks = fitting.PeakFinderFitter.getConvergedPeaks(self, peaks)
        return self.peak_fitter.rescaleZ(converged_peaks)


