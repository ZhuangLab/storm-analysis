#!/usr/bin/python
#
# Cubic spline peak finder.
#
# Hazen 01/14
#

import pickle
import numpy

import sa_library.fitting as fitting
import sa_library.ia_utilities_c as util_c

import cubic_fit_c as cubic_fit_c


#
# Spliner peak finding.
#
class SplinerPeakFinder(fitting.PeakFinder):
    
    def findPeaks(self, image, peaks):
        # Set peak finding cutoff based on current background and threshold.
        self.cutoff = self.background + self.cur_threshold

        # Find the peaks using the standard peak finder.
        return fitting.PeakFinder.findPeaks(self, image, peaks)

#
# Spliner peak fitting.
#
class SplinerPeakFitter(fitting.PeakFitter):

    def __init__(self, fitting_function, parameters):
        fitting.PeakFitter.__init__(self, fitting_function, parameters)

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
            self.sfitter = cubic_fit_c.CSpline2DFit(self.spline, self.coeff, False)
        else:
            self.spline_type = "3D"
            self.sfitter = cubic_fit_c.CSpline3DFit(self.spline, self.coeff, False)

        # Save the coefficients for faster start up.
        if save_coeff:
            psf_data["coeff"] = self.sfitter.getCoeff()
            pickle.dump(psf_data, open(parameters.spline, "w"))

    def fitPeaks(self, peaks):

        # Fit to update peak locations.
        self.sfitter.doFit(peaks)
        fit_peaks = self.sfitter.getGoodPeaks(min_height = 0.9*self.threshold)

        # Remove peaks that are too close to each other & refit.
        fit_peaks = util_c.removeClosePeaks(fit_peaks, self.sigma, self.neighborhood)

        # Redo the fit for the remaining peaks.
        self.sfitter.doFit(fit_peaks)
        fit_peaks = self.sfitter.getGoodPeaks(min_height = 0.9*self.threshold)
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
        self.peak_finder = SplinerPeakFinder(parameters, self.margin)
        self.peak_fitter = SplinerPeakFitter(False, parameters)

    def analyzeImage(self, new_image, save_residual = False, verbose = False):
        return fitting.PeakFinderFitter.analyzeImage(self, new_image, save_residual = True)

    def getConvergedPeaks(self, peaks):
        converged_peaks = fitting.PeakFinderFitter.getConvergedPeaks(self, peaks)
        return self.peak_fitter.rescaleZ(converged_peaks)


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
