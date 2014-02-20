#!/usr/bin/python
#
# sCMOS peak finder.
#
# Hazen 01/14
#

import numpy

import sa_library.fitting as fitting
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

    def __init__(self, parameters, margin):
        fitting.PeakFinder.__init__(self, parameters, margin)
        
        # Create image smoother object.
        [lg_offset, lg_variance, lg_gain] = loadSCMOSData(parameters.camera_calibration,
                                                          margin)
        self.smoother = scmos_utilities_c.Smoother(lg_offset, lg_variance, lg_gain)

    def findPeaks(self, image, peaks):
        # Calculate convolved image.
        smooth_image = self.smoother.smoothImage(image, self.sigma)

        # Set peak finding cutoff. Since finding is performed on a convolved
        # image it is not necessary to include a background offset.
        self.cutoff = self.cur_threshold

        # Find the peaks using the standard peak finder.
        return fitting.PeakFinder.findPeaks(self, smooth_image, peaks)


#
# sCMOS peak fitting.
#
class SCMOSPeakFitter(fitting.PeakFitter):

    def __init__(self, fitting_function, parameters, margin):
        fitting.PeakFitter.__init__(self, fitting_function, parameters)

        # Create image regularizer object & calibration term for peak fitting.
        [lg_offset, lg_variance, lg_gain] = loadSCMOSData(parameters.camera_calibration,
                                                          margin)
        self.scmos_cal = lg_variance/(lg_gain*lg_gain)
        self.regularizer = scmos_utilities_c.Regularizer(lg_offset, lg_variance, lg_gain)

    def fitPeaks(self, peaks):
        [fit_peaks, residual] = fitting.PeakFitter.fitPeaks(self, peaks)
        residual = self.regularizer.deregularizeImage(residual)
        
        return [fit_peaks, residual]
        
    def newImage(self, new_image):
        self.image = self.regularizer.regularizeImage(new_image)


#
# Base class to encapsulate sCMOS peak finding and fitting.
#
class SCMOSFinderFitter(fitting.PeakFinderFitter):

    def __init__(self, parameters):
        fitting.PeakFinderFitter.__init__(self, parameters)
        self.peak_finder = SCMOSPeakFinder(parameters, self.margin)


class SCMOS2DFixed(SCMOSFinderFitter):

    def __init__(self, parameters):
        SCMOSFinderFitter.__init__(self, parameters)
        self.peak_fitter = SCMOSPeakFitter(multi_c.fitMultiGaussian2DFixed,
                                           parameters,
                                           self.margin)


class SCMOS2D(SCMOSFinderFitter):

    def __init__(self, parameters):
        SCMOSFinderFitter.__init__(self, parameters)
        self.peak_fitter = SCMOSPeakFitter(multi_c.fitMultiGaussian2D,
                                           parameters,
                                           self.margin)


class SCMOS3D(SCMOSFinderFitter):

    def __init__(self, parameters):
        SCMOSFinderFitter.__init__(self, parameters)
        self.peak_fitter = SCMOSPeakFitter(multi_c.fitMultiGaussian3D,
                                           parameters,
                                           self.margin)


class SCMOSZ(SCMOSFinderFitter):

    def __init__(self, parameters):
        SCMOSFinderFitter.__init__(self, parameters)
        self.peak_fitter = SCMOSPeakFitter(multi_c.fitMultiGaussianZ,
                                           parameters,
                                           self.margin)


#
# Return the appropriate type of fitter.
#
def initFindAndFit(parameters):
    options = [["2dfixed", SCMOS2DFixed],
               ["2d", SCMOS2D],
               ["3d", SCMOS3D],
               ["Z", SCMOSZ]]
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
