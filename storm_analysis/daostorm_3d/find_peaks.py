#!/usr/bin/env python
"""

Finds "all" the peaks in an image.

Hazen 01/14
"""

import numpy

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.matched_filter_c as matchedFilterC
import storm_analysis.sa_library.dao_fit_c as daoFitC
import storm_analysis.simulator.draw_gaussians_c as dg


class Daostorm2DFixedFitter(fitting.PeakFitter):
    def __init__(self, **kwds):
        super(Daostorm2DFixedFitter, self).__init__(**kwds)
        self.mfitter = daoFitC.MultiFitter2DFixed(self.scmos_cal, self.wx_params, self.wy_params, self.min_z, self.max_z)


class Daostorm2DFitter(fitting.PeakFitter):

    def __init__(self, **kwds):
        super(Daostorm2DFitter, self).__init__(**kwds)
        self.mfitter = daoFitC.MultiFitter2D(self.scmos_cal, self.wx_params, self.wy_params, self.min_z, self.max_z)
        

class Daostorm3DFitter(fitting.PeakFitter):

    def __init__(self, **kwds):
        super(Daostorm2DFixedFitter, self).__init__(**kwds)
        self.mfitter = daoFitC.MultiFitter3D(self.scmos_cal, self.wx_params, self.wy_params, self.min_z, self.max_z)
        
    
class DaostormZFitter(fitting.PeakFitter):
    
    def __init__(self, **kwds):
        super(Daostorm2DFixedFitter, self).__init__(**kwds)
        self.mfitter = daoFitC.MultiFitterZ(self.scmos_cal, self.wx_params, self.wy_params, self.min_z, self.max_z)
        

def initFindAndFit(parameters):
    """
    Return the appropriate type of finder and fitter.
    """
    finder = fitting.PeakFinder(parameters = parameters)

    # Initialize fitter.
    fitters = {'2dfixed' : Daostorm2DFixedFitter,
               '2d' : Daostorm2DFitter,
               '3d' : Daostorm3DFitter,
               'Z' : DaostormZFitter}
    fitter = fitters[parameters.getAttr("model")](parameters = parameters)
    
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
