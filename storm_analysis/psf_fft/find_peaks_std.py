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


class PSFFFTPeakFinder(fitting.PeakFinderArbitraryPSF):
    """
    PSF FFT peak finding.
    """
    def __init__(self, parameters = None, psf_fn = None, **kwds):
        kwds["parameters"] = parameters
        super(PSFFFTPeakFinder, self).__init__(**kwds)
        self.psf_object = psf_fn

        # Update margin based on the PSF FFT size.
        old_margin = self.margin
        self.margin = int((self.psf_object.getSize() + 1)/2 + 2)

        self.fg_mfilter_zval = parameters.getAttr("z_value", [0.0])
        for zval in self.fg_mfilter_zval:
            self.z_values.append(self.psf_object.getScaledZ(zval))

        # Load peak locations if specified.
        #
        # Note: This is not in the base class because different fitters use different scales
        #       for the Z parameter.
        #
        if parameters.hasAttr("peak_locations"):
            [self.peak_locations, is_text] = fitting.getPeakLocations(parameters.getAttr("peak_locations"),
                                                                      self.margin,
                                                                      parameters.getAttr("pixel_size"),
                                                                      self.sigma)

            zc_index = utilC.getZCenterIndex()
            # Set initial z value (for text files).
            if is_text:
                self.peak_locations[:,zc_index] = self.z_value[0]

            # Convert z value to PSF FFT units (Insight3 localization files).
            else:
                for i in range(self.peak_locations.shape[0]):
                    self.peak_locations[i,zc_index] = self.psf_object.getScaledZ(self.peak_locations[i,zc_index])


class PSFFFTPeakFitter(fitting.PeakFitter):
    """
    PSF FFT peak fitting.
    """
    def __init__(self, **kwds):
        super(PSFFFTPeakFitter, self).__init__(**kwds)

        # Update refitting neighborhood parameter.
        self.neighborhood = int(0.25 * self.mfitter.getSize()) + 1
        
    # Convert from spline z units to real z units.
    def rescaleZ(self, peaks):
        return self.mfitter.rescaleZ(peaks)
        

class PSFFFTPeakFinderFitter(fitting.PeakFinderFitter):
    """
    Class for spline based peak finding and fitting.
    """
    def getConvergedPeaks(self, peaks):
        converged_peaks = super(PSFFFTPeakFinderFitter, self).getConvergedPeaks(peaks)
        return self.peak_fitter.rescaleZ(converged_peaks)


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
    finder = PSFFFTPeakFinder(parameters = parameters,
                              psf_fn = psf_fn)

    # Create fftFitC.CFFTFit object.
    mfitter = initFitter(finder, parameters, psf_fn)
    
    # Create peak fitter.
    fitter = PSFFFTPeakFitter(mfitter = mfitter,
                              parameters = parameters)

    return PSFFFTPeakFinderFitter(peak_finder = finder,
                                  peak_fitter = fitter)

