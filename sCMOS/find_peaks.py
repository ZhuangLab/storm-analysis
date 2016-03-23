#!/usr/bin/python
#
# sCMOS peak finder.
#
# Hazen 01/14
#

import numpy

import sa_library.fitting as fitting
import sa_library.ia_utilities_c as util_c
import sa_library.multi_fit_c as multi_c

import scmos_utilities_c


# Load sCMOS data
def loadSCMOSData(calibration_filename, margin):

    # Additional sCMOS specific data & objects.
    # Load camera calibrations & create smoother and regularizer.
    [offset, variance, gain] = numpy.load(calibration_filename)

    # Pad out sCMOS arrays.
    lg_offset = fitting.padArray(offset, margin)
    lg_variance = fitting.padArray(variance, margin)
    lg_gain = fitting.padArray(gain, margin)

    return [lg_offset, lg_variance, lg_gain]


#
# sCMOS peak finding.
#
class SCMOSPeakFinder(fitting.PeakFinder):

    def __init__(self, parameters):
        fitting.PeakFinder.__init__(self, parameters)
        
        # Create image smoother object.
        [lg_offset, lg_variance, lg_gain] = loadSCMOSData(parameters.camera_calibration,
                                                          fitting.PeakFinderFitter.margin)
        self.smoother = scmos_utilities_c.Smoother(lg_offset, lg_variance, lg_gain)

    def peakFinder(self, image):

        # Calculate convolved image, this is where the background subtraction happens.
        smooth_image = self.smoother.smoothImage(image, self.sigma)

        # Mask the image so that peaks are only found in the AOI.
        masked_image = smooth_image * self.peak_mask
        
        # Identify local maxima in the masked image.
        [new_peaks, self.taken] = util_c.findLocalMaxima(masked_image,
                                                         self.taken,
                                                         self.cur_threshold,
                                                         self.find_max_radius,
                                                         self.margin)
        return new_peaks

    def subtractBackground(self, image):

        # Estimate the background from the current (residual) image.
        self.background = self.backgroundEstimator(image)

        #
        # Just return the image as peakFinder() will do the background subtraction by
        # convolution and we don't want to make things complicated by also doing that here.
        #
        return image
    
    
#
# sCMOS peak fitting.
#
class SCMOSPeakFitter(fitting.PeakFitter):

    def __init__(self, parameters):
        fitting.PeakFitter.__init__(self, parameters)

        # Create image regularizer object & calibration term for peak fitting.
        [lg_offset, lg_variance, lg_gain] = loadSCMOSData(parameters.camera_calibration,
                                                          fitting.PeakFinderFitter.margin)
        self.scmos_cal = lg_variance/(lg_gain*lg_gain)
        self.regularizer = scmos_utilities_c.Regularizer(lg_offset, lg_variance, lg_gain)

    def fitPeaks(self, peaks):
        [fit_peaks, residual] = fitting.PeakFitter.fitPeaks(self, peaks)
        residual = self.regularizer.deregularizeImage(residual)
        
        return [fit_peaks, residual]
        
    def newImage(self, new_image):
        self.image = self.regularizer.regularizeImage(new_image)


class SCMOS2DFixedFitter(SCMOSPeakFitter):

    def peakFitter(self, peaks):
        return multi_c.fitMultiGaussian2DFixed(self.image, peaks, self.scmos_cal)


class SCMOS2DFitter(SCMOSPeakFitter):

    def peakFitter(self, peaks):
        return multi_c.fitMultiGaussian2D(self.image, peaks, self.scmos_cal)


class SCMOS3DFitter(SCMOSPeakFitter):

    def peakFitter(self, peaks):
        return multi_c.fitMultiGaussian3D(self.image, peaks, self.scmos_cal)

    
class SCMOSZFitter(SCMOSPeakFitter):

    def peakFitter(self, peaks):
        return multi_c.fitMultiGaussianZ(self.image, peaks, self.scmos_cal)
        

#
# Base class to encapsulate sCMOS peak finding and fitting.
#
class sCMOSFinderFitter(fitting.PeakFinderFitter):

    def __init__(self, parameters, peak_finder, peak_fitter):
        fitting.PeakFinderFitter.__init__(self, parameters)
        self.peak_finder = peak_finder
        self.peak_fitter = peak_fitter
        

#
# Return the appropriate type of finder and fitter.
#
def initFindAndFit(parameters):
    fitters = {'2dfixed' : SCMOS2DFixedFitter,
               '2d' : SCMOS2DFitter,
               '3d' : SCMOS3DFitter,
               'Z' : SCMOSZFitter}
    return sCMOSFinderFitter(parameters,
                             SCMOSPeakFinder(parameters),
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
