#!/usr/bin/env python
"""
Initialize for multi-plane peak finding and fitting with arbitrary PSFs.

Hazen 01/18
"""
import os
import pickle
import numpy
import tifffile

import storm_analysis.multi_plane.fitting_mp as fittingMp
import storm_analysis.multi_plane.mp_fit_arb_c as mpFitArbC
import storm_analysis.multi_plane.mp_utilities as mpUtil

import storm_analysis.sa_library.analysis_io as analysisIO

import storm_analysis.psf_fft.psf_fn as psfFn
import storm_analysis.pupilfn.pupil_fn as pupilFn
import storm_analysis.spliner.spline_to_psf as splineToPSF

    
def initFitter(margin, parameters, psf_objects, rqes, variances):
    """
    Create and return a mpFitArbC.MPXXFit object.
    """
    # This fitter only supports 'MLE'
    assert (parameters.getAttr("fit_error_model") == 'MLE'), "Only MLE fitting is supported."
    
    assert(len(psf_objects) == len(variances))
    #
    # FIXME: Not sure having two z ranges, one from the spline
    #        and one in the parameters is a good idea.
    #
    # Create the fitter object which will do the actual fitting. Unless specified
    # the fit for each channel is forced to have the same height.
    #
    if isinstance(psf_objects[0], psfFn.PSFFn):
        mfitter = mpFitArbC.MPPSFFnFit(independent_heights = parameters.getAttr("independent_heights", 0),
                                       psf_objects = psf_objects)

    elif isinstance(psf_objects[0], pupilFn.PupilFunction):
        mfitter = mpFitArbC.MPPupilFnFit(independent_heights = parameters.getAttr("independent_heights", 0),
                                         psf_objects = psf_objects)

    elif isinstance(psf_objects[0], splineToPSF.SplineToPSF3D):
        mfitter = mpFitArbC.MPSplineFit(independent_heights = parameters.getAttr("independent_heights", 0),
                                        psf_objects = psf_objects)

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

    # Load channel Cramer-Rao weights if available.
    #
    weights = None
    if parameters.hasAttr("weights"):
        with open(parameters.getAttr("weights"), 'rb') as fp:
            weights = pickle.load(fp)
    mfitter.setWeights(weights)
    
    return mfitter


def initPSFObjects(parameters):
    """
    Create and return the PSF objects (spline, pupil function or psf FFT).
    """
    psf_objects = []

    # Try PSF FFT.
    #
    if (len(mpUtil.getPSFFFTAttrs(parameters)) > 0):
        
        # Create PSF FFT PSF objects.
        for psf_fft_attr in mpUtil.getPSFFFTAttrs(parameters):
            psf_objects.append(psfFn.PSFFn(psf_filename = parameters.getAttr(psf_fft_attr)))

        # All the PSF FFT objects have to have the same Z range.
        for i in range(1, len(psf_objects)):
            assert (psf_objects[0].getZMin() == psf_objects[i].getZMin())
            assert (psf_objects[0].getZMax() == psf_objects[i].getZMax())
            
    # Try pupil functions.
    #
    elif (len(mpUtil.getPupilFnAttrs(parameters)) > 0):

        # Get fitting Z range (this is in microns).
        [min_z, max_z] = parameters.getZRange()
        
        # Create pupil function PSF objects.
        for pupil_fn_attr in mpUtil.getPupilFnAttrs(parameters):
            psf_objects.append(pupilFn.PupilFunction(pf_filename = parameters.getAttr(pupil_fn_attr),
                                                     zmin = min_z * 1.0e+3,
                                                     zmax = max_z * 1.0e+3))

    # Try splines.
    #
    elif (len(mpUtil.getSplineAttrs(parameters)) > 0):

        # Create Spline PSF objects.
        for spline_attr in mpUtil.getSplineAttrs(parameters):
            psf_objects.append(splineToPSF.SplineToPSF3D(spline_file = parameters.getAttr(spline_attr)))

        # All the splines have to have the same Z range.
        for i in range(1, len(psf_objects)):
            assert (psf_objects[0].getZMin() == psf_objects[i].getZMin())
            assert (psf_objects[0].getZMax() == psf_objects[i].getZMax())

    else:
        raise Exception("No PSF objects found.")

    return psf_objects


def initFindAndFit(parameters):
    """
    Create and return a fittingMp.MPFinderFitter object.
    """
    # Create PSF objects.
    psf_objects = initPSFObjects(parameters)
    
    # Create peak finder.
    finder = fittingMp.MPPeakFinderArb(parameters = parameters,
                                       psf_objects = psf_objects)

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
                         psf_objects,
                         rqes,
                         variances)

    # Create peak fitter.
    #
    fitter = fittingMp.MPPeakFitterArb(mfitter = mfitter,
                                       parameters = parameters)

    # Specify which properties we want (for each channel) from the
    # analysis. Note that some of these may be duplicates of each
    # other, for example if the heights are not independent.
    #
    properties = ["background", "error", "height", "iterations", "significance", "sum", "x", "y", "z"]

    return fittingMp.MPFinderFitter(peak_finder = finder,
                                    peak_fitter = fitter,
                                    properties = properties)
