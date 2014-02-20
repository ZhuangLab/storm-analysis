#!/usr/bin/python
#
# Finds "all" the peaks in an image.
#
# Hazen 01/14
#

import numpy

import sa_library.fitting as fitting
import sa_library.multi_fit_c as multi_c

#
# Daostorm peak finding.
#
class DaostormPeakFinder(fitting.PeakFinder):
    
    def findPeaks(self, image, peaks):
        # Set peak finding cutoff based on current background and threshold.
        self.cutoff = self.background + self.cur_threshold

        # Find the peaks using the standard peak finder.
        return fitting.PeakFinder.findPeaks(self, image, peaks)


#
# Base class to encapsulate 3d-daostorm peak finding and fitting.
#
class DaostormFinderFitter(fitting.PeakFinderFitter):

    def __init__(self, parameters):
        fitting.PeakFinderFitter.__init__(self, parameters)
        self.peak_finder = DaostormPeakFinder(parameters, self.margin)
        

class Daostorm2DFixed(DaostormFinderFitter):
    
    def __init__(self, parameters):
        DaostormFinderFitter.__init__(self, parameters)
        self.peak_fitter = fitting.PeakFitter(multi_c.fitMultiGaussian2DFixed, parameters)


class Daostorm2D(DaostormFinderFitter):
    
    def __init__(self, parameters):
        DaostormFinderFitter.__init__(self, parameters)
        self.peak_fitter = fitting.PeakFitter(multi_c.fitMultiGaussian2D, parameters)


class Daostorm3D(DaostormFinderFitter):
    
    def __init__(self, parameters):
        DaostormFinderFitter.__init__(self, parameters)
        self.peak_fitter = fitting.PeakFitter(multi_c.fitMultiGaussian3D, parameters)


class DaostormZ(DaostormFinderFitter):

    def __init__(self, parameters):
        DaostormFinderFitter.__init__(self, parameters)
        self.peak_fitter = fitting.PeakFitter(multi_c.fitMultiGaussianZ, parameters)


#
# Return the appropriate type of fitter.
#
def initFindAndFit(parameters):
    options = [["2dfixed", Daostorm2DFixed],
               ["2d", Daostorm2D],
               ["3d", Daostorm3D],
               ["Z", DaostormZ]]
    return fitting.initFindAndFit(parameters, options)

#
# The MIT License
#
# Copyright (c) 2014 Zhuang Lab, Harvard University
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
