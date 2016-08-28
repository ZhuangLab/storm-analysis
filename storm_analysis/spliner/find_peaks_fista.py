#!/usr/bin/python
#
# FISTA deconvolution based peak finder.
# Cubic spline based fitting.
#
# FIXME: This doesn't support 2D spline fitting.
#
#
# Hazen 01/16
#

import pickle
import numpy

import storm_analysis.fista.fista_decon as fistaDecon
import storm_analysis.rolling_ball_bgr.rolling_ball as rollingBall
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.wavelet_bgr.wavelet_bgr as waveletBGR

import cubic_fit_c as cubicFitC
import spline_to_psf as splineToPSF

#
# FISTA peak finding.
#
class SplinerPeakFinder(object):

    def __init__(self, parameters):
        self.fista_iterations = parameters.fista_iterations
        self.fista_lambda = parameters.fista_lambda
        self.fista_number_z = parameters.fista_number_z
        self.fista_threshold = parameters.fista_threshold
        self.fista_timestep = parameters.fista_timestep
        self.fista_upsample = parameters.fista_upsample
        self.spline_file = parameters.spline

        self.rball = None
        self.wbgr = None

        # Update margin based on the spline size.
        s_to_psf = splineToPSF.SplineToPSF(parameters.spline)
        self.margin = (s_to_psf.getSize() + 1)/4 + 2
        
        if hasattr(parameters, "rb_radius"):
            self.rball = rollingBall.RollingBall(parameters.rb_radius,
                                                 parameters.rb_sigma)
        else:
            self.wbgr_iterations = parameters.wbgr_iterations
            self.wbgr_threshold = parameters.wbgr_threshold
            self.wbgr_wavelet_level = parameters.wbgr_wavelet_level
            self.wbgr = waveletBGR.WaveletBGR()
            
        self.fdecon = None
    
    def findPeaks(self):
        
        # Run the FISTA deconvolution.
        self.fdecon.decon(self.fista_iterations, self.fista_lambda)

        # Get the peaks from the deconvolved image.
        peaks = self.fdecon.getPeaks(self.fista_threshold, self.margin)

        return peaks

    def newImage(self, image, bg_estimate):

        # Create FISTA deconvolver if it doesn't exist.
        if self.fdecon is None:
            self.fdecon = fistaDecon.FISTADecon(image.shape,
                                                self.spline_file,
                                                self.fista_number_z,
                                                self.fista_timestep,
                                                self.fista_upsample)
        
            print("Margin is", self.margin)

        # Use provided background estimate.
        #
        # FIXME: Should we run background subtraction on the estimated background?
        #
        if bg_estimate is not None:
            background = bg_estimate
            
        # Estimate background.
        else:
            if self.rball is not None:
                # Use rolling ball approach.
                background = self.rball.estimateBG(image)
            
            else:
                # Use wavelet approach.
                background = self.wbgr.estimateBG(image,
                                                  self.wbgr_iterations,
                                                  self.wbgr_threshold,
                                                  self.wbgr_wavelet_level)

        # Configure FISTA solver with the new image & estimated background.
        self.fdecon.newImage(image, background)

        
#
# Spliner peak fitting.
#
class SplinerPeakFitter(object):

    def __init__(self, parameters):
        self.fit_threshold = parameters.fit_threshold
        self.fit_sigma = parameters.fit_sigma

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
            print("2D spline fitting is not supported.")
            exit()
            #self.spline_type = "2D"
            #self.sfitter = cubic_fit_c.CSpline2DFit(self.spline, self.coeff, False)
        else:
            self.spline_type = "3D"
            self.sfitter = cubicFitC.CSpline3DFit(self.spline, self.coeff, False)

        # Save the coefficients for faster start up.
        if save_coeff:
            psf_data["coeff"] = self.sfitter.getCoeff()
            pickle.dump(psf_data, open(parameters.spline, "w"))

        # Calculate refitting neighborhood parameter.
        self.fit_neighborhood = int(0.25 * self.sfitter.getSize()) + 1
        
    def fitPeaks(self, peaks):

        # Adjust to z starting position.
        z_index = utilC.getZCenterIndex()
        peaks[:,z_index] = peaks[:,z_index] * float(self.sfitter.getSize())

        if False:
            print("Before fitting")
            for i in range(5):
                print(" ", peaks[i,0], peaks[i,1], peaks[i,3], peaks[i,5], peaks[i,6], peaks[i,7])
            print("")

        # Fit to update peak locations.
        self.sfitter.doFit(peaks)
        fit_peaks = self.sfitter.getGoodPeaks(min_height = 0.9 * self.fit_threshold)

        # Remove peaks that are too close to each other & refit.
        fit_peaks = utilC.removeClosePeaks(fit_peaks, self.fit_sigma, self.fit_neighborhood)

        # Redo the fit for the remaining peaks.
        self.sfitter.doFit(fit_peaks)
        fit_peaks = self.sfitter.getGoodPeaks(min_height = 0.9*self.fit_threshold)
        residual = self.sfitter.getResidual()

        if False:
            print("After fitting")
            for i in range(5):
                print(" ", fit_peaks[i,0], fit_peaks[i,1], fit_peaks[i,3], fit_peaks[i,5], fit_peaks[i,6], fit_peaks[i,7])
            print("")
        
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
class SplinerFinderFitter(object):

    def __init__(self, parameters):
        self.peak_finder = SplinerPeakFinder(parameters)
        self.peak_fitter = SplinerPeakFitter(parameters)

    #
    # FIXME:
    #   bg_estimate handling has not been tested.
    #
    def analyzeImage(self, new_image, bg_estimate = None, save_residual = False, verbose = False):
        
        image = fitting.padArray(new_image, self.peak_finder.margin)
        if bg_estimate is not None:
            bg_estimate = fitting.padArray(bg_estimate, self.peak_finder.margin)
            
        self.peak_finder.newImage(image, bg_estimate)
        self.peak_fitter.newImage(image)
        
        #
        # This is a lot simpler than 3D-DAOSTORM as we only do one pass,
        # hopefully the compressed sensing (FISTA) deconvolution finds all the
        # peaks and then we do a single pass of fitting.
        #
        if True:
            peaks = self.peak_finder.findPeaks()
            [fit_peaks, residual] = self.peak_fitter.fitPeaks(peaks)

        #
        # This is for testing if just using FISTA followed by the center
        # of mass calculation is basically as good as also doing the
        # additional MLE spline fitting step.
        #
        # The short answer is that it appears that it is not. It about
        # 1.3x worse in XY and about 4x worse in Z.
        #
        else:
            fit_peaks = self.peak_finder.findPeaks()

            # Adjust z scale.
            z_index = utilC.getZCenterIndex()
            z_size = (self.peak_fitter.spline.shape[2] - 1.0)
            status_index = utilC.getStatusIndex()
            fit_peaks[:,z_index] = z_size*fit_peaks[:,z_index]
            
            # Mark as converged.
            fit_peaks[:,status_index] = 1.0
            
            residual = None

        #
        # Subtract margin so that peaks are in the right
        # place with respect to the original image.
        #
        fit_peaks[:,utilC.getXCenterIndex()] -= float(self.peak_finder.margin)
        fit_peaks[:,utilC.getYCenterIndex()] -= float(self.peak_finder.margin)

        return [fit_peaks, residual]

    def cleanUp(self):
        pass

    def getConvergedPeaks(self, peaks, verbose = False):
        if (peaks.shape[0] > 0):
            status_index = utilC.getStatusIndex()
            mask = (peaks[:,status_index] == 1.0)  # 0.0 = running, 1.0 = converged.
            if verbose:
                print(" ", numpy.sum(mask), "converged out of", peaks.shape[0])
            return self.peak_fitter.rescaleZ(peaks[mask,:])
        else:
            return peaks
    

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
