#!/usr/bin/env python
"""
PSF FFT peak finder / fitter.

Hazen 10/17
"""

import pickle
import numpy

import tifffile

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
    # Variance is in units of ADU*ADU.
    # Gain is ADU/photo-electron.
    #
    variance = None
    if parameters.hasAttr("camera_calibration"):
        [offset, variance, gain] = numpy.load(parameters.getAttr("camera_calibration"))
        variance = variance/(gain*gain)

        # Set variance in the peak finder, this method also pads the
        # variance to the correct size.
        variance = finder.setVariance(variance)

    # Create C fitter object.
    return fftFitC.CFFTFit(psf_fn = psf_fn,
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

    return fitting.PeakFinderFitterArbitraryPSF(peak_finder = finder,
                                                peak_fitter = fitter)

