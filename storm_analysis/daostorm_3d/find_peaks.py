#!/usr/bin/env python
"""
Finds "all" the peaks in an image.

Hazen 01/14
"""
import numpy

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.dao_fit_c as daoFitC

def initFindAndFit(parameters):
    """
    Return the appropriate type of finder and fitter.
    """
    fmodel = parameters.getAttr("model")
    emodel = parameters.getAttr("fit_error_model")
    
    # Create peak finder.
    finder = fitting.PeakFinderGaussian(parameters = parameters)

    # Initialize Z fitting parameters.
    wx_params = None
    wy_params = None
    min_z = None
    max_z = None
    if (parameters.getAttr("model", "na") == "Z"):
        [wx_params, wy_params] = parameters.getWidthParams(for_mu_Zfit = True)
        [min_z, max_z] = parameters.getZRange()

    # Check for camera calibration (this function is also used by sCMOS analysis).
    rqe = None
    variance = None
    if parameters.hasAttr("camera_calibration"):
        
        # Load variance, scale by gain.
        #
        # Offset is in units of ADU.
        # Variance is in units of ADU*ADU.
        # Gain is ADU/photo-electron.
        # RQE is dimensionless, it should be around 1.0.
        #
        [offset, variance, gain, rqe] = analysisIO.loadCMOSCalibration(parameters.getAttr("camera_calibration"))
        variance = variance/(gain*gain)

        # Set variance in the peak finder, this method also pads the
        # variance to the correct size.
        variance = finder.setVariance(variance)

        # Pad relative quantum efficiency array to the correct size.
        rqe = finder.padArray(rqe)

    # Create C fitter object.
    mfitter = None
    kwds = {'roi_size' : finder.getROISize(),
            'rqe' : rqe,
            'scmos_cal' : variance,
            'wx_params' : wx_params,
            'wy_params' : wy_params,
            'min_z' : min_z,
            'max_z' : max_z}
    
    if (fmodel == '2dfixed'):
        if (emodel == 'MLE'):
            mfitter = daoFitC.MultiFitter2DFixed(**kwds)
        elif (emodel == 'ALS'):
            mfitter = daoFitC.MultiFitter2DFixedALS(**kwds)
        elif (emodel == 'LS'):
            mfitter = daoFitC.MultiFitter2DFixedLS(**kwds)
        elif (emodel == 'DWLS'):
            mfitter = daoFitC.MultiFitter2DFixedDWLS(**kwds)
        elif (emodel == 'FWLS'):
            mfitter = daoFitC.MultiFitter2DFixedFWLS(**kwds)

        # These models are for tests.
        elif (emodel == 'MLE_NC'):
            mfitter = daoFitC.MultiFitter2DFixedNC(**kwds)
        elif (emodel == 'DWLS_NC'):
            mfitter = daoFitC.MultiFitter2DFixedDWLSNC(**kwds)
        
    elif (fmodel == '2d'):
        if parameters.hasAttr("sigma_range"):
            kwds['sigma_range'] = parameters.getAttr("sigma_range")
        else:
            sigma = parameters.getAttr("sigma")
            kwds['sigma_range'] = [0.5 * sigma, 5.0 * sigma]
            
        if (emodel == 'MLE'):
            mfitter = daoFitC.MultiFitter2D(**kwds)
        elif (emodel == 'ALS'):
            mfitter = daoFitC.MultiFitter2DALS(**kwds)
        elif (emodel == 'LS'):
            mfitter = daoFitC.MultiFitter2DLS(**kwds)
        elif (emodel == 'DWLS'):
            mfitter = daoFitC.MultiFitter2DDWLS(**kwds)
        elif (emodel == 'FWLS'):
            mfitter = daoFitC.MultiFitter2DFWLS(**kwds)

    elif (fmodel == '3d'):
        if parameters.hasAttr("sigma_range"):
            kwds['sigma_range'] = parameters.getAttr("sigma_range")
        else:
            sigma = parameters.getAttr("sigma")
            kwds['sigma_range'] = [0.5 * sigma, 5.0 * sigma]
            
        if (emodel == 'MLE'):
            mfitter = daoFitC.MultiFitter3D(**kwds)

    elif (fmodel == 'Z'):
        if (emodel == 'MLE'):
            mfitter = daoFitC.MultiFitterZ(**kwds)

    if mfitter is None:
        raise Exception("The requested fitting model and/or error model is not available. '" + fmodel + "' '" + emodel + "'")

    # Create peak fitter.
    fitter = fitting.PeakFitter(mfitter = mfitter,
                                parameters = parameters)

    # Specify which properties we want from the analysis.
    properties = ["background", "error", "height", "iterations", "significance", "sum", "x", "y"]
    if (fmodel == "2dfixed") or (fmodel == "2d"):
        properties.append("xsigma")
    elif (fmodel == "3d"):
        properties.append("xsigma")
        properties.append("ysigma")
    elif (fmodel == "Z"):
        properties.append("xsigma")
        properties.append("ysigma")
        properties.append("z")

    return fitting.PeakFinderFitter(peak_finder = finder,
                                    peak_fitter = fitter,
                                    properties = properties)

#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
