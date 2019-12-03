#!/usr/bin/env python
"""
CS deconvolution based peak finder.
Cubic spline based fitting.

FIXME: This doesn't support 2D spline fitting.

Hazen 11/19
"""
import pickle
import numpy
import warnings

import storm_analysis

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.parameters as params

import storm_analysis.rolling_ball_bgr.rolling_ball as rollingBall
import storm_analysis.wavelet_bgr.wavelet_bgr as waveletBGR

import storm_analysis.admm.admm_decon as admmDecon
import storm_analysis.densestorm.densestorm_decon as densestormDecon
import storm_analysis.fista.fista_decon as fistaDecon

import storm_analysis.spliner.cubic_fit_c as cubicFitC
import storm_analysis.spliner.find_peaks_std as findPeaksStd
import storm_analysis.spliner.spline_to_psf as splineToPSF


class FindPeaksDeconException(storm_analysis.SAException):
    pass


class CSDeconPeakFinder(object):
    """
    Base class for compressed sensing deconvolution peak finding.
    """
    def __init__(self, bg_estimator = None, parameters = None, psf_object = None, **kwds):
        super(CSDeconPeakFinder, self).__init__(**kwds)

        self.bg_estimator = bg_estimator
        self.decon_object = None
        self.psf_object = psf_object
        
        self.iterations = parameters.getAttr("iterations")
        if (self.iterations != 1):
            warnings.warn("Peak finding iterations is > 1 for deconvolution!")

        # Update margin based on the psf object size.
        self.margin = self.psf_object.getMargin()
            
    def cleanUp(self):
        if self.decon_object is not None:
            self.decon_object.cleanup()
            self.decon_object = None

    def estimateBackground(self, fit_peaks_image, bg_estimate):
        """
        This gets called once per cycle of peak finding and fitting. At each
        cycle fit_peaks_image is the current best estimate of the foreground
        only portion of the image.
        """
        
        # Use provided background estimate.
        if bg_estimate is not None:
            self.background = bg_estimate
            
        # Estimate background.
        else:
            image = self.image - fit_peaks_image
            self.background = self.bg_estimator.estimateBG(image)

        # Pass new background estimate to the CS decon solver.
        self.decon_object.newBackground(self.background)
        
        return self.background

    def findPeaks(self, fit_peaks_image):
        """
        Find the peaks in the image.
        """
        # Run the deconvolution.
        self.deconvolve()

        # Get the peaks from the deconvolved image.
        peaks = self.getPeaks()

        # Convert z values from nanometers to PSF units.
        peaks["z"] = self.psf_object.getScaledZ(peaks["z"])
        
        return [peaks, "finder", True]

    def getDWLSError(self):
        """
        Return CS solvers current DWLS error.
        """
        return self.decon_object.getDWLSError()
    
    def getL1Error(self):
        """
        Return CS solvers current l1 error.
        """
        return self.decon_object.getL1Error()

    def getL2Error(self):
        """
        Return CS solvers current l2 error.
        """
        return self.decon_object.getL2Error()

    def getXVector(self):
        """
        Return x vector (the deconvolved image) from the CS solver.
        """
        return self.decon_object.getXVector()
        
    def newImage(self, image):
        """
        Setup to segment a new image. This gets called once per image.
        """
        self.image = numpy.copy(image)

        # Create CS solver if it doesn't exist. We do this here
        # because we need to know the dimensions of the image
        # in order to do this.
        if self.decon_object is None:
            self.deconInit(image)

        # Pass new image to the CS solver.
        self.decon_object.newImage(image)

    def setVariance(self, camera_variance):
        """
        Just return the camera_variance array properly re-sized.
        """
        return fitting.padArray(camera_variance, self.margin)


class ADMMPeakFinder(CSDeconPeakFinder):
    """
    ADMM deconvolution peak finding.
    """
    def __init__(self, parameters = None, **kwds):
        super(ADMMPeakFinder, self).__init__(parameters = parameters, **kwds)

        self.admm_iterations = parameters.getAttr("admm_iterations")
        self.admm_lambda = parameters.getAttr("admm_lambda")
        self.admm_number_z = parameters.getAttr("admm_number_z")
        self.admm_rho = parameters.getAttr("admm_rho")
        self.admm_threshold = parameters.getAttr("admm_threshold")

    def deconInit(self, image):
        self.decon_object = admmDecon.ADMMDecon(image.shape,
                                                self.psf_object,
                                                self.admm_number_z,
                                                self.admm_rho)

    def deconvolve(self):
        self.decon_object.decon(self.admm_iterations, self.admm_lambda)

    def getPeaks(self):
        return self.decon_object.getPeaks(self.admm_threshold, self.margin)

    
class DenseSTORMPeakFinder(CSDeconPeakFinder):
    """
    3denseSTORM deconvolution peak finding.
    """
    def __init__(self, parameters = None, **kwds):
        super(DenseSTORMPeakFinder, self).__init__(parameters = parameters, **kwds)

        self.ds_beta = parameters.getAttr("ds3_beta")
        self.ds_eta = parameters.getAttr("ds3_eta")
        self.ds_iterations = parameters.getAttr("ds3_iterations")
        self.ds_micro = parameters.getAttr("ds3_micro")
        self.ds_number_z = parameters.getAttr("ds3_number_z")
        self.ds_threshold = parameters.getAttr("ds3_threshold")

    def deconInit(self, image):
        self.decon_object = densestormDecon.DenseSTORMDecon(image.shape,
                                                            self.psf_object,
                                                            self.ds_number_z,
                                                            self.ds_beta,
                                                            self.ds_eta,
                                                            self.ds_micro)

    def deconvolve(self):
        self.decon_object.decon(self.ds_iterations)

    def getPeaks(self):
        return self.decon_object.getPeaks(self.ds_threshold, self.margin)

    
class FISTAPeakFinder(CSDeconPeakFinder):
    """
    FISTA deconvolution peak finding.
    """
    def __init__(self, parameters = None, **kwds):
        super(FISTAPeakFinder, self).__init__(parameters = parameters, **kwds)

        self.fista_iterations = parameters.getAttr("fista_iterations")
        self.fista_lambda = parameters.getAttr("fista_lambda")
        self.fista_number_z = parameters.getAttr("fista_number_z")
        self.fista_threshold = parameters.getAttr("fista_threshold")
        self.fista_timestep = parameters.getAttr("fista_timestep")

    def deconInit(self, image):
        self.decon_object = fistaDecon.FISTADecon(image.shape,
                                                  self.psf_object,
                                                  self.fista_number_z,
                                                  self.fista_timestep)

    def deconvolve(self):
        self.decon_object.decon(self.fista_iterations, self.fista_lambda)

    def getPeaks(self):
        return self.decon_object.getPeaks(self.fista_threshold, self.margin)


def initBGEstimator(parameters, settings_name):
    """
    Return requested background estimator object.

    Note: As a side effect this will add additional values to parameters.
    """
    estimator_name = parameters.getAttr("background_estimator")

    bg_estimator = None
    if (estimator_name == "RollingBall"):
        parameters.update(params.ParametersRollingBall().initFromFile(settings_name, warnings = False))
        bg_estimator = rollingBall.RollingBall(parameters.getAttr("rb_radius"),
                                               parameters.getAttr("rb_sigma"))

    elif (estimator_name == "Wavelet"):
        parameters.update(params.ParametersWaveletBGR().initFromFile(settings_name, warnings = False))
        bg_estimator = waveletBGR.WaveletBGRStdAna(iterations = parameters.getAttr("wbgr_iterations"),
                                                   threshold = parameters.getAttr("wbgr_threshold"),
                                                   wavelet_level = parameters.getAttr("wbgr_wavelet_level"))

    else:
        raise FindPeaksDeconException("Unknown background estimator " + estimator_name)

    return bg_estimator


def initDeconFinder(parameters, settings_name, psf_object):
    """
    Return a decon peak finding object.

    Note: As a side effect this will add additional values to parameters.
    """
    finder = None
    
    if (parameters.getAttr("decon_method") == "ADMM"):
        parameters.update(params.ParametersADMM().initFromFile(settings_name, warnings = False))
        bg_estimator = initBGEstimator(parameters, settings_name)
        finder = ADMMPeakFinder(bg_estimator = bg_estimator,
                                parameters = parameters,
                                psf_object = psf_object)

    elif (parameters.getAttr("decon_method") == "3denseSTORM"):
        parameters.update(params.Parameters3denseSTORM().initFromFile(settings_name, warnings = False))
        bg_estimator = initBGEstimator(parameters, settings_name)
        finder = DenseSTORMPeakFinder(bg_estimator = bg_estimator,
                                      parameters = parameters,
                                      psf_object = psf_object)

    elif (parameters.getAttr("decon_method") == "FISTA"):
        parameters.update(params.ParametersFISTA().initFromFile(settings_name, warnings = False))
        bg_estimator = initBGEstimator(parameters, settings_name)
        finder = FISTAPeakFinder(bg_estimator = bg_estimator,
                                 parameters = parameters,
                                 psf_object = psf_object)

    else:
        raise FindPeaksDeconException("Unknown deconvolution method " + parameters.getAttr("decon_method"))

    return finder


def initFindAndFit(parameters, settings_name):
    """
    Initialize and return a PeakFinderFitter object.
    """
    # Create spline object.
    spline_fn = splineToPSF.loadSpline(parameters.getAttr("spline"))

    # Create peak finder.
    finder = initDeconFinder(parameters, settings_name, spline_fn)
    
    # Create cubicFitC.CSplineFit object.
    mfitter = findPeaksStd.initFitter(finder, parameters, spline_fn)
    
    # Create peak fitter.
    fitter = fitting.PeakFitterArbitraryPSF(mfitter = mfitter,
                                            parameters = parameters)

    # Specify which properties we want from the analysis.
    properties = ["background", "error", "height", "iterations", "significance", "sum", "x", "y", "z"]
    
    return fitting.PeakFinderFitter(peak_finder = finder,
                                    peak_fitter = fitter,
                                    properties = properties)

