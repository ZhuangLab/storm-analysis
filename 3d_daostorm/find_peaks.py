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
# 3d-daostorm peak finding.
#
class DaostormPeakFinder(fitting.PeakFinder):
    pass


#
# 3d-daostorm peak fitting.
#
class Daostorm2DFixedFitter(fitting.PeakFitter):

    def peakFitter(self, peaks):
        return multi_c.fitMultiGaussian2DFixed(self.image, peaks, None)


class Daostorm2DFitter(fitting.PeakFitter):

    def peakFitter(self, peaks):
        return multi_c.fitMultiGaussian2D(self.image, peaks, None)


class Daostorm3DFitter(fitting.PeakFitter):

    def peakFitter(self, peaks):
        return multi_c.fitMultiGaussian3D(self.image, peaks, None)

    
class DaostormZFitter(fitting.PeakFitter):

    def peakFitter(self, peaks):
        return multi_c.fitMultiGaussianZ(self.image, peaks, None)


#
# Base class to encapsulate 3d-daostorm peak finding and fitting.
#
class DaostormFinderFitter(fitting.PeakFinderFitter):

    def __init__(self, parameters, peak_finder, peak_fitter):
        fitting.PeakFinderFitter.__init__(self, parameters)
        self.peak_finder = peak_finder
        self.peak_fitter = peak_fitter
        

#
# Return the appropriate type of finder and fitter.
#
def initFindAndFit(parameters):
    fitters = {'2dfixed' : Daostorm2DFixedFitter,
               '2d' : Daostorm2DFitter,
               '3d' : Daostorm3DFitter,
               'Z' : DaostormZFitter}
    return DaostormFinderFitter(parameters,
                                DaostormPeakFinder(parameters),
                                fitters[parameters.model](parameters))

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
