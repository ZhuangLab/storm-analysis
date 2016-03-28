#!/usr/bin/python
#
# Finds "all" the peaks in an image.
#
# Hazen 01/14
#

import numpy

import sa_library.fitting as fitting
import sa_library.ia_utilities_c as utilC
import sa_library.matched_filter_c as matchedFilterC
import sa_library.multi_fit_c as multiC
import simulator.drawgaussians as dg

#
# 3D-DAOSTORM peak finding (for low SNR data).
#
class DaostormPeakFinder(fitting.PeakFinder):

    def __init__(self, parameters):
        fitting.PeakFinder.__init__(self, parameters)        
        self.filter_sigma = parameters.filter_sigma
        self.mfilter = None

    def newImage(self, new_image):
        fitting.PeakFinder.newImage(self, new_image)

        # If does not already exist, create a gaussian filter object.
        if self.mfilter is None:
            psf = dg.drawGaussiansXY(new_image.shape,
                                     numpy.array([0.5*new_image.shape[0]]),
                                     numpy.array([0.5*new_image.shape[1]]),
                                     sigma = self.filter_sigma)
            psf = psf/numpy.sum(psf)
            self.mfilter = matchedFilterC.MatchedFilter(psf)

    def peakFinder(self, image):
        
        # Smooth image with gaussian filter.
        smooth_image = self.mfilter.convolve(image)
        
        # Mask the image so that peaks are only found in the AOI.
        masked_image = smooth_image * self.peak_mask
        
        # Identify local maxima in the masked image.
        [new_peaks, self.taken] = utilC.findLocalMaxima(masked_image,
                                                        self.taken,
                                                        self.cur_threshold,
                                                        self.find_max_radius,
                                                        self.margin)
        return new_peaks


#
# 3D-DAOSTORM peak fitting.
#
class Daostorm2DFixedFitter(fitting.PeakFitter):

    def peakFitter(self, peaks):
        return multiC.fitMultiGaussian2DFixed(self.image, peaks, None)


class Daostorm2DFitter(fitting.PeakFitter):

    def peakFitter(self, peaks):
        return multiC.fitMultiGaussian2D(self.image, peaks, None)


class Daostorm3DFitter(fitting.PeakFitter):

    def peakFitter(self, peaks):
        return multiC.fitMultiGaussian3D(self.image, peaks, None)

    
class DaostormZFitter(fitting.PeakFitter):

    def peakFitter(self, peaks):
        return multiC.fitMultiGaussianZ(self.image, peaks, None)


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

    # Initialize finder.
    if hasattr(parameters, "filter_sigma") and (parameters.filter_sigma > 0.0):
        finder = DaostormPeakFinder(parameters)
    else:
        finder = fitting.PeakFinder(parameters)

    # Initialize fitter.
    fitters = {'2dfixed' : Daostorm2DFixedFitter,
               '2d' : Daostorm2DFitter,
               '3d' : Daostorm3DFitter,
               'Z' : DaostormZFitter}
    fitter = fitters[parameters.model](parameters)
    
    return DaostormFinderFitter(parameters, finder, fitter)

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
