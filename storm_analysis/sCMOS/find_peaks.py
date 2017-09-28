#!/usr/bin/env python
"""
sCMOS peak finder.

Hazen 01/14
"""

import numpy

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.dao_fit_c as daoFitC


class SCMOSPeakFitter(fitting.PeakFitter):
    """
    sCMOS peak fitting.
    """
    def __init__(self, camera_variance = None, **kwds):
        super(SCMOSPeakFitter, self).__init__(**kwds)

        self.scmos_cal = camera_variance

class SCMOS2DFixedFitter(SCMOSPeakFitter):
    
    def __init__(self, **kwds):
        super(SCMOS2DFixedFitter, self).__init__(**kwds)
        self.mfitter = daoFitC.MultiFitter2DFixed(self.scmos_cal, self.wx_params, self.wy_params, self.min_z, self.max_z)


class SCMOS2DFitter(SCMOSPeakFitter):

    def __init__(self, **kwds):
        super(SCMOS2DFitter, self).__init__(**kwds)
        self.mfitter = daoFitC.MultiFitter2D(self.scmos_cal, self.wx_params, self.wy_params, self.min_z, self.max_z)


class SCMOS3DFitter(SCMOSPeakFitter):

    def __init__(self, **kwds):
        super(SCMOS3DFitter, self).__init__(**kwds)
        self.mfitter = daoFitC.MultiFitter3D(self.scmos_cal, self.wx_params, self.wy_params, self.min_z, self.max_z)

    
class SCMOSZFitter(SCMOSPeakFitter):

    def __init__(self, **kwds):
        super(SCMOSZFitter, self).__init__(**kwds)
        self.mfitter = daoFitC.MultiFitterZ(self.scmos_cal, self.wx_params, self.wy_params, self.min_z, self.max_z)


def initFindAndFit(parameters):
    """
    Return the appropriate type of finder and fitter.
    """
    # Create peak finder.
    finder = fitting.PeakFinder(parameters = parameters)

    # Load variance, scale by gain.
    #
    # Variance is in units of ADU*ADU.
    # Gain is ADU/photo-electron.
    #
    [offset, variance, gain] = numpy.load(parameters.getAttr("camera_calibration"))
    variance = variance/(gain*gain)

    # Set variance in the peak finder, this method also pads the
    # variance to the correct size.
    variance = finder.setVariance(variance)

    fitters = {'2dfixed' : SCMOS2DFixedFitter,
               '2d' : SCMOS2DFitter,
               '3d' : SCMOS3DFitter,
               'Z' : SCMOSZFitter}
    fitter = fitters[parameters.getAttr("model")](parameters = parameters,
                                                  camera_variance = variance)

    return fitting.PeakFinderFitter(peak_finder = finder,
                                    peak_fitter = fitter)


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
