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
    def __init__(self, parameters = None, psf_object = None, **kwds):
        super(SplinerFISTAPeakFinder, self).__init__(**kwds)

        self.psf_object = psf_object
        
        self.fista_iterations = parameters.getAttr("fista_iterations")
        self.fista_lambda = parameters.getAttr("fista_lambda")
        self.fista_number_z = parameters.getAttr("fista_number_z")
        self.fista_threshold = parameters.getAttr("fista_threshold")
        self.fista_timestep = parameters.getAttr("fista_timestep")

        # Only perform a single cycle of peak finding and fitting.
        self.iterations = 1

        self.rball = None
        self.wbgr = None

        # Update margin based on the spline size.
        self.margin = self.psf_object.getMargin()

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
        self.fdecon.cleanup()
    
    def findPeaks(self, fit_peaks_image):
        """
        Find the peaks in the image.
        """        
        # Run the FISTA deconvolution.
        self.fdecon.decon(self.fista_iterations, self.fista_lambda)

        # Get the peaks from the deconvolved image.
        peaks = self.fdecon.getPeaks(self.fista_threshold, self.margin)

        return [peaks, "finder", True]

    def newImage(self, image):
        """
        Setup to segment a new image, mostly this involves background estimation.
        """
        # Create FISTA deconvolver if it doesn't exist.
        if self.fdecon is None:
            self.fdecon = fistaDecon.FISTADecon(image.shape,
                                                self.psf_object,
                                                self.fista_number_z,
                                                self.fista_timestep)

        # Pass new image to the FISTA solver.
        self.fdecon.newImage(image)

    def setVariance(self, camera_variance):
        """
        Just return the camera_variance array properly re-sized.
        """
        return fitting.padArray(camera_variance, self.margin)

    def subtractBackground(self, image, bg_estimate):
        
        # Use provided background estimate.
        if bg_estimate is not None:
            self.background = bg_estimate
            
        # Estimate background.
        else:
            if self.rball is not None:
                # Use rolling ball approach.
                self.background = self.rball.estimateBG(image)
                
            else:
                # Use wavelet approach.
                self.background = self.wbgr.estimateBG(image,
                                                       self.wbgr_iterations,
                                                       self.wbgr_threshold,
                                                       self.wbgr_wavelet_level)

        # Pass new background estimate to the FISTA solver.
        self.fdecon.newBackground(self.background)
        
        return self.background

    
def initFindAndFit(parameters):
    """
    Initialize and return a SplinerFISTAFinderFitter object.
    """
    # Create spline object.
    spline_fn = splineToPSF.loadSpline(parameters.getAttr("spline"))    
    
    # Create peak finder.
    finder = SplinerFISTAPeakFinder(parameters = parameters,
                                    psf_object = spline_fn)

    # Create cubicFitC.CSplineFit object.
    mfitter = findPeaksStd.initFitter(finder, parameters, spline_fn)
    
    # Create peak fitter.
    fitter = fitting.PeakFitterArbitraryPSF(mfitter = mfitter,
                                            parameters = parameters)

    # Specify which properties we want from the analysis.
    properties = ["background", "error", "height", "x", "y", "z"]
    
    return fitting.PeakFinderFitter(peak_finder = finder,
                                    peak_fitter = fitter,
                                    properties = properties)
