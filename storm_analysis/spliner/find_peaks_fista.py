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

import storm_analysis.spliner.cubic_fit_c as cubicFitC
import storm_analysis.spliner.spline_to_psf as splineToPSF


class FindPeaksFistaException(Exception):
    pass

#
# Spliner FISTA peak finding.
#
class SplinerFISTAPeakFinder(object):

    def __init__(self, parameters):
        self.fista_iterations = parameters.getAttr("fista_iterations")
        self.fista_lambda = parameters.getAttr("fista_lambda")
        self.fista_number_z = parameters.getAttr("fista_number_z")
        self.fista_threshold = parameters.getAttr("fista_threshold")
        self.fista_timestep = parameters.getAttr("fista_timestep")
        self.fista_upsample = parameters.getAttr("fista_upsample")
        self.spline_file = parameters.getAttr("spline")

        self.rball = None
        self.wbgr = None

        # Load spline to get size.
        self.spline_file = parameters.getAttr("spline")
        s_to_psf = splineToPSF.loadSpline(self.spline_file)

        # Update margin based on the spline size.
        self.margin = int((s_to_psf.getSize() + 1)/4 + 2)
        
        if parameters.hasAttr("rb_radius"):
            self.rball = rollingBall.RollingBall(parameters.getAttr("rb_radius"),
                                                 parameters.getAttr("rb_sigma"))
        else:
            self.wbgr_iterations = parameters.getAttr("wbgr_iterations")
            self.wbgr_threshold = parameters.getAttr("wbgr_threshold")
            self.wbgr_wavelet_level = parameters.getAttr("wbgr_wavelet_level")
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
class SplinerFISTAPeakFitter(object):

    def __init__(self, parameters):
        self.threshold = parameters.getAttr("threshold")
        self.sigma = parameters.getAttr("sigma")

        # Load spline and create the appropriate type of spline fitter.
        psf_data = pickle.load(open(parameters.getAttr("spline"), 'rb'))
        if (psf_data["type"] == "3D"):
            self.zmin = psf_data["zmin"]/1000.0
            self.zmax = psf_data["zmax"]/1000.0
        self.spline = psf_data["spline"]
        self.coeff = psf_data["coeff"]
            
        if (len(self.spline.shape) == 2):
            self.spline_type = "2D"
            self.mfitter = cubicFitC.CSpline2DFit(self.spline, self.coeff, None)
        else:
            self.spline_type = "3D"
            self.mfitter = cubicFitC.CSpline3DFit(self.spline, self.coeff, None)

        # Calculate refitting neighborhood parameter.
        self.fit_neighborhood = int(0.25 * self.mfitter.getSplineSize()) + 1

    def fitPeaks(self, peaks):

        # Adjust to z starting position.
        z_index = utilC.getZCenterIndex()
        peaks[:,z_index] = peaks[:,z_index] * float(self.mfitter.getSplineSize())

        if False:
            print("Before fitting")
            for i in range(5):
                print(" ", peaks[i,0], peaks[i,1], peaks[i,3], peaks[i,5], peaks[i,6], peaks[i,7])
            print("")

        # Fit to update peak locations.
        fit_peaks = self.mfitter.doFit(peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks, min_height = 0.9 * self.threshold)

        # Remove peaks that are too close to each other & refit.
        fit_peaks = utilC.removeClosePeaks(fit_peaks, self.sigma, self.fit_neighborhood)

        # Redo the fit for the remaining peaks.
        fit_peaks = self.mfitter.doFit(fit_peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks, min_height = 0.9*self.threshold)
        residual = self.mfitter.getResidual()

        if False:
            print("After fitting")
            for i in range(5):
                print(" ", fit_peaks[i,0], fit_peaks[i,1], fit_peaks[i,3], fit_peaks[i,5], fit_peaks[i,6], fit_peaks[i,7])
            print("")
        
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
# Spline fitting using FISTA for peak finding.
#
class SplinerFISTAFinderFitter(object):

    def __init__(self, parameters):
        self.peak_finder = SplinerFISTAPeakFinder(parameters)
        self.peak_fitter = SplinerFISTAPeakFitter(parameters)

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
    

