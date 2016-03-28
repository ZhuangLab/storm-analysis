#!/usr/bin/python
#
# Cubic spline peak finder.
#
# Hazen 03/16
#

import pickle
import numpy

import tifffile

import sa_library.fitting as fitting
import sa_library.ia_utilities_c as utilC
import sa_library.matched_filter_c as matchedFilterC

import cubic_fit_c as cubicFitC
import spline_to_psf as splineToPSF


#
# Spliner peak finding.
#
class SplinerPeakFinder(fitting.PeakFinder):

    def __init__(self, parameters):
        fitting.PeakFinder.__init__(self, parameters)        
        self.s_to_psf = splineToPSF.SplineToPSF(parameters.spline)
        self.mfilter = None

    def newImage(self, new_image):
        fitting.PeakFinder.newImage(self, new_image)

        # If does not already exist, create a gaussian filter object
        # from the best fit spline to the PSF.
        if self.mfilter is None:
            psf = self.s_to_psf.getPSF(0, shape = new_image.shape)
            self.mfilter = matchedFilterC.MatchedFilter(psf)

            # Save a picture of the PSF for debugging purposes.
            if 0:
                temp = 10000.0 * psf
                tifffile.imsave("psf.tif", temp.astype(numpy.uint16))
            
    def peakFinder(self, image):
        
        # Smooth image with gaussian filter.
        smooth_image = self.mfilter.convolve(image)
        
        # Mask the image so that peaks are only found in the AOI.
        masked_image = smooth_image * self.peak_mask
        
        # Identify local maxima in the masked image.
        [new_peaks, self.taken] = utilC.findLocalMaxima(masked_image,
                                                        self.taken,
                                                        self.cur_threshold,
                                                        self.find_max_radius,
                                                        self.margin)
        
        return new_peaks


#
# Spliner peak fitting.
#
class SplinerPeakFitter(fitting.PeakFitter):

    def __init__(self, parameters):
        fitting.PeakFitter.__init__(self, parameters)

        # Load spline and create the appropriate type of spline fitter.
        psf_data = pickle.load(open(parameters.spline))
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
            self.sfitter = cubicFitC.CSpline2DFit(self.spline, self.coeff, False)
        else:
            self.spline_type = "3D"
            self.sfitter = cubicFitC.CSpline3DFit(self.spline, self.coeff, False)

        # Save the coefficients for faster start up.
        if save_coeff:
            psf_data["coeff"] = self.sfitter.getCoeff()
            pickle.dump(psf_data, open(parameters.spline, "w"))

    def fitPeaks(self, peaks):

        # Fit to update peak locations.
        self.sfitter.doFit(peaks)
        fit_peaks = self.sfitter.getGoodPeaks(min_height = 0.9 * self.threshold,
                                              verbose = False)

        # Remove peaks that are too close to each other & refit.
        fit_peaks = utilC.removeClosePeaks(fit_peaks, self.sigma, self.neighborhood)

        # Redo the fit for the remaining peaks.
        self.sfitter.doFit(fit_peaks)
        fit_peaks = self.sfitter.getGoodPeaks(min_height = 0.9*self.threshold,
                                              verbose = False)
        residual = self.sfitter.getResidual()

        return [fit_peaks, residual]

    def newImage(self, image):
        self.sfitter.newImage(image)

    # Convert from spline z units to real z units.
    def rescaleZ(self, peaks):
        if (self.spline_type == "3D"):
            return self.sfitter.rescaleZ(peaks, self.zmin, self.zmax)
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

    def analyzeImage(self, new_image, save_residual = False, verbose = False):
        return fitting.PeakFinderFitter.analyzeImage(self, new_image, save_residual = True)

    def getConvergedPeaks(self, peaks):
        converged_peaks = fitting.PeakFinderFitter.getConvergedPeaks(self, peaks)
        return self.peak_fitter.rescaleZ(converged_peaks)


#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
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
