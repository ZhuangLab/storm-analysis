#!/usr/bin/env python
"""
PSF FFT peak finder / fitter.

Hazen 10/17
"""

import pickle
import numpy

import tifffile

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC

import storm_analysis.psf_fft.fft_fit_c as fftFitC
import storm_analysis.psf_fft.psf_fn as psfFn


def initFitter(finder, parameters, psf_fn):
    """
    Initialize and return a fftFitC.CFFTFit object.
    """
    # Load variance, scale by gain.
    #
    # Offset is in units of ADU.
    # Variance is in units of ADU*ADU.
    # Gain is ADU/photo-electron.
    # RQE is dimensionless, it should be around 1.0.
    #

    # This fitter only supports 'MLE'
    assert (parameters.getAttr("fit_error_model") == 'MLE'), "Only MLE fitting is supported."
    
    rqe = None
    variance = None
    if parameters.hasAttr("camera_calibration"):
        [offset, variance, gain, rqe] = analysisIO.loadCMOSCalibration(parameters.getAttr("camera_calibration"))
        variance = variance/(gain*gain)

        # Set variance in the peak finder, this method also pads the
        # variance to the correct size.
        variance = finder.setVariance(variance)
        
        # Pad relative quantum efficiency array to the correct size.
        rqe = finder.padArray(rqe)

    # Create C fitter object.
    return fftFitC.CFFTFit(psf_fn = psf_fn,
                           rqe = rqe,
                           scmos_cal = variance)


def initFindAndFit(parameters):
    """
    Initialize and return a PSFFFTFinderFitter object.
    """
    # Create psf fft function object.
    psf_fn = psfFn.PSFFn(psf_filename = parameters.getAttr("psf"))

    # Check that the PSF FFT and camera pixel sizes agree.
    diff = abs(parameters.getAttr("pixel_size") - psf_fn.getPixelSize()*1.0e3)
    assert (diff < 1.0e-6), "Wrong pixel size, incorrect PSF data?"
    
    # Create peak finder.
    finder = fitting.PeakFinderArbitraryPSF(parameters = parameters,
                                            psf_object = psf_fn)

    # Create fftFitC.CFFTFit object.
    mfitter = initFitter(finder, parameters, psf_fn)
    
    # Create peak fitter.
    fitter = fitting.PeakFitterArbitraryPSF(mfitter = mfitter,
                                            parameters = parameters)
    
    # Specify which properties we want from the analysis.
    properties = ["background", "error", "height", "iterations", "significance", "sum", "x", "y", "z"]

    return fitting.PeakFinderFitter(peak_finder = finder,
                                    peak_fitter = fitter,
                                    properties = properties)

