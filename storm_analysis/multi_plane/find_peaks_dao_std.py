#!/usr/bin/env python
"""
Initialize for multi-plane peak finding and fitting with Gaussian PSFs.

Hazen 01/18
"""
import os
import pickle
import numpy
import tifffile

import storm_analysis.multi_plane.fitting_mp as fittingMp
import storm_analysis.multi_plane.mp_utilities as mpUtil
import storm_analysis.sa_library.analysis_io as analysisIO

import storm_analysis.multi_plane.mp_fit_dao_c as mpFitDaoC


def initFitter(margin, parameters, roi_size, rqes, variances):
    """
    Create and return a mpFitDaoC.MPFitDao2D object.
    """
    # This fitter only supports 'MLE'
    assert (parameters.getAttr("fit_error_model") == 'MLE'), "Only MLE fitting is supported."
    
    sigma = parameters.getAttr("sigma")
    mfitter = mpFitDaoC.MPFitDao2D(n_channels = len(variances),
                                   roi_size = roi_size,
                                   sigma_range = [0.5 * sigma, 3.0 * sigma])

    # Initialize fitters for each channel.
    #
    for i in range(len(variances)):
        mfitter.initializeChannel(rqes[i], variances[i], i)

    # Load mappings.
    #
    if parameters.hasAttr("mapping"):
        if os.path.exists(parameters.getAttr("mapping")):
            mapping_filename = parameters.getAttr("mapping")
        else:
            raise Exception("Mapping file", parameters.getAttr("mapping"), "does not exist.")

        mfitter.setMapping(*mpUtil.loadMappings(mapping_filename, margin))

    # Initialize weights.
    #
    mfitter.setWeights()
    
    return mfitter


def initFindAndFit(parameters):
    """
    Create and return a fittingMp.MPFinderFitter object.
    """
    # Load sCMOS calibration data.
    #
    # Note: Gain is expected to be in units of ADU per photo-electron.
    #
    rqes = []
    variances = []
    for calib_name in mpUtil.getCalibrationAttrs(parameters):
        [offset, variance, gain, rqe] = analysisIO.loadCMOSCalibration(parameters.getAttr(calib_name))
        variances.append(variance/(gain*gain))
        rqes.append(rqe)
    
    # Create peak finder.
    finder = fittingMp.MPPeakFinderDao(parameters = parameters,
                                       n_channels = len(variances))

    # Set variance in the peak finder. This method also pads the variance
    # to the correct size and performs additional initializations.
    #
    variances = finder.setVariances(variances)

    # Pad the rqes to the correct size.
    rqes = list(map(lambda x: finder.padArray(x), rqes))

    # Create mpFitC.MPFit object.
    #
    mfitter = initFitter(finder.margin,
                         parameters,
                         finder.roi_size,
                         rqes,
                         variances)

    # Create peak fitter.
    #
    fitter = fittingMp.MPPeakFitterDao(mfitter = mfitter,
                                       parameters = parameters)

    # Specify which properties we want (for each channel) from the
    # analysis. Note that some of these may be duplicates of each
    # other, for example if the heights are not independent.
    #
    properties = ["background", "error", "height", "iterations", "significance", "sum", "x", "xsigma", "y"]

    return fittingMp.MPFinderFitter(peak_finder = finder,
                                    peak_fitter = fitter,
                                    properties = properties)
