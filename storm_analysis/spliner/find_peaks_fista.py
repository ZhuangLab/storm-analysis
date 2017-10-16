#!/usr/bin/env python
"""
FISTA deconvolution based peak finder.
Cubic spline based fitting.

FIXME: This doesn't support 2D spline fitting.

Hazen 01/16
"""

import pickle
import numpy

import storm_analysis.fista.fista_decon as fistaDecon
import storm_analysis.rolling_ball_bgr.rolling_ball as rollingBall
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.wavelet_bgr.wavelet_bgr as waveletBGR

import storm_analysis.spliner.cubic_fit_c as cubicFitC
import storm_analysis.spliner.find_peaks_std as findPeaksStd
import storm_analysis.spliner.spline_to_psf as splineToPSF


class FindPeaksFistaException(Exception):
    pass


class SplinerFISTAPeakFinder(object):
    """
    Spliner FISTA peak finding.
    """
    def __init__(self, parameters = None, **kwds):
        super(SplinerFISTAPeakFinder, self).__init__(**kwds)
        
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

    def cleanUp(self):
        pass
    
    def findPeaks(self):
        """
        Find the peaks in the image.
        """        
        # Run the FISTA deconvolution.
        self.fdecon.decon(self.fista_iterations, self.fista_lambda)

        # Get the peaks from the deconvolved image.
        peaks = self.fdecon.getPeaks(self.fista_threshold, self.margin)

        return peaks

    def newImage(self, image, bg_estimate):
        """
        Setup to segment a new image, mostly this involves background estimation.
        """
        # Create FISTA deconvolver if it doesn't exist.
        if self.fdecon is None:
            self.fdecon = fistaDecon.FISTADecon(image.shape,
                                                self.spline_file,
                                                self.fista_number_z,
                                                self.fista_timestep,
                                                self.fista_upsample)
        
            print("Margin is", self.margin)

        # Use provided background estimate.
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

    def setVariance(self, camera_variance):
        """
        Just return the camera_variance array properly re-sized.
        """
        return fitting.padArray(camera_variance, self.margin)


class SplinerFISTAPeakFitter(findPeaksStd.SplinerPeakFitter):
    """
    Spliner peak fitting.
    """
    def fitPeaks(self, peaks):

        # Adjust to z starting position.
        z_index = utilC.getZCenterIndex()
        peaks[:,z_index] = peaks[:,z_index] * float(self.mfitter.getSize())

        if False:
            print("Before fitting")
            for i in range(5):
                print(" ", peaks[i,0], peaks[i,1], peaks[i,3], peaks[i,5], peaks[i,6], peaks[i,7])
            print("")

        # Fit to update peak locations.
        fit_peaks = self.mfitter.doFit(peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks, 0.0)

        # Remove peaks that are too close to each other & refit.
        fit_peaks = utilC.removeClosePeaks(fit_peaks, self.sigma, self.neighborhood)

        # Redo the fit for the remaining peaks.
        fit_peaks = self.mfitter.doFit(fit_peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks, 0.0)
        fit_peaks_image = self.mfitter.getResidual()

        if False:
            print("After fitting")
            for i in range(5):
                print(" ", fit_peaks[i,0], fit_peaks[i,1], fit_peaks[i,3], fit_peaks[i,5], fit_peaks[i,6], fit_peaks[i,7])
            print("")
        
        return [fit_peaks, fit_peaks_image]
        

class SplinerFISTAFinderFitter(findPeaksStd.SplinerFinderFitter):
    """
    Spline fitting using FISTA for peak finding.
    """
    #
    # FIXME: bg_estimate handling has not been tested.
    #
    def analyzeImage(self, movie_reader, save_residual = False, verbose = False):

        # Load image (in photo-electrons).
        [image, fit_peaks_image] = self.loadImage(movie_reader)
        
        # Load background estimate (in photo-electrons).
        bg_estimate = self.loadBackgroundEstimate(movie_reader)
            
        self.peak_finder.newImage(image, bg_estimate)
        self.peak_fitter.newImage(image)
        
        #
        # This is a lot simpler than 3D-DAOSTORM as we only do one pass,
        # hopefully the compressed sensing (FISTA) deconvolution finds all the
        # peaks and then we do a single pass of fitting.
        #
        if True:
            peaks = self.peak_finder.findPeaks()
            [fit_peaks, fit_peaks_image] = self.peak_fitter.fitPeaks(peaks)

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

        return [fit_peaks, fit_peaks_image]

    
def initFindAndFit(parameters):
    """
    Initialize and return a SplinerFISTAFinderFitter object.
    """
    # Create peak finder.
    finder = SplinerFISTAPeakFinder(parameters = parameters)

    # Create cubicFitC.CSplineFit object.
    mfitter = findPeaksStd.initFitter(finder, parameters)
    
    # Create peak fitter.
    fitter = SplinerFISTAPeakFitter(mfitter = mfitter,
                                    parameters = parameters)

    return SplinerFISTAFinderFitter(peak_finder = finder,
                                    peak_fitter = fitter)    
