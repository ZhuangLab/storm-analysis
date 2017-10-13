#!/usr/bin/env python
"""
Pupil function peak fitting.

Hazen 10/17
"""

import pickle
import numpy

import tifffile

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC

import storm_analysis.pupilfn.pupil_fit_c as pupilFitC
import storm_analysis.pupilfn.pupil_fn as pupilFn


class PupilFnPeakFinder(fitting.PeakFinderArbitraryPSF):
    """
    Pupil function peak finding.
    """
    def __init__(self, parameters = None, pupil_fn = None, **kwds):
        kwds["parameters"] = parameters
        super(PupilFnPeakFinder, self).__init__(**kwds)
        
        # Load the pupil function.
        self.psf_object = pupil_fn

        # Update margin based on the pupil function size (this
        # is always an even number).
        old_margin = self.margin
        self.margin = int(self.psf_object.getSize()/2 + 2)

        self.fg_mfilter_zval = parameters.getAttr("z_value", [0.0])
        self.z_values = self.fg_mfilter_zval

        if parameters.hasAttr("peak_locations"):

            # Correct for any difference in the margins.
            self.peak_locations[:,utilC.getXCenterIndex()] += self.margin - old_margin
            self.peak_locations[:,utilC.getYCenterIndex()] += self.margin - old_margin

            # Provide the "correct" starting z value.
            self.peak_locations[:,utilC.getZCenterIndex()] = self.z_values[0]


class PupilFnPeakFitter(fitting.PeakFitter):
    """
    Pupil function peak fitting.
    """
    def __init__(self, **kwds):
        super(PupilFnPeakFitter, self).__init__(**kwds)

        # Update refitting neighborhood parameter.
        self.neighborhood = int(0.5 * self.mfitter.getSize()) + 1


def initFitter(finder, parameters, pupil_fn):
    """
    Initialize and return a pupilFitC.CPupilFit object.
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

    # Get fitting Z range.
    [min_z, max_z] = parameters.getZRange()

    # Create C fitter object.
    return pupilFitC.CPupilFit(pupil_fn = pupil_fn,
                               scmos_cal = variance,
                               min_z = min_z,
                               max_z = max_z)


def initFindAndFit(parameters):
    """
    Initialize and return a SplinerFinderFitter object.
    """
    # Create pupil function object.
    pupil_fn = pupilFn.PupilFunction(pf_filename = parameters.getAttr("pupil_function"))

    # Check that the PF and camera pixel sizes agree.
    diff = abs(parameters.getAttr("pixel_size") - pupil_fn.getPixelSize()*1.0e3)
    assert (diff < 1.0e-6), "Incorrect pupil function?"
    
    # Create peak finder.
    finder = PupilFnPeakFinder(parameters = parameters,
                               pupil_fn = pupil_fn)

    # Create cubicFitC.CSplineFit object.
    mfitter = initFitter(finder, parameters, pupil_fn)
    
    # Create peak fitter.
    fitter = PupilFnPeakFitter(mfitter = mfitter,
                               parameters = parameters)

    return fitting.PeakFinderFitter(peak_finder = finder,
                                    peak_fitter = fitter)
