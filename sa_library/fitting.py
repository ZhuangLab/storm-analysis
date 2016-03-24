#!/usr/bin/python
#
## @file
#
# This contains some of what is common to all of the peak finding algorithms.
#
# FitData: A base class for storing all the parameters associated with
#    image analysis.
#
# FindAndFit: A base class for peak finding and fitting.
#
# Hazen 01/14
#

import numpy
import scipy
import scipy.ndimage

import daxwriter
import ia_utilities_c as util_c
import multi_fit_c as multi_c
import parameters as params


#
# Functions.
#

## estimateBackground
#
# A simple background estimator.
#
# @param image A 2D numpy array containing the image.
# @param size (Optional) The size of the filter to use, default is 8.
#
# @return A 2D numpy array containing the image background estimate.
#
def estimateBackground(image, size = 8):
    background = scipy.ndimage.filters.gaussian_filter(image, (size, size))
    return background

## initZParams
#
# Initialize parameters for Z fitting.
#
# @param wx_params Array of wx parameters.
# @param wy_params Array of wy parameters.
# @param min_z Minimum allowable z value.
# @param nax_z Maximum allowable z value.
#
def initZParams(wx_params, wy_params, min_z, max_z):
    multi_c.initZParams(wx_params, wy_params, min_z, max_z)

## padArray
#
# Pads out an array to a large size.
#
# @param ori_array A 2D numpy array.
# @param pad_size The number of elements to add to each of the "sides" of the array.
#
# @return The padded 2D numpy array.
#
def padArray(ori_array, pad_size):
    [x_size, y_size] = ori_array.shape
    lg_array = numpy.ones((x_size+2*pad_size,y_size+2*pad_size))
    lg_array[pad_size:(x_size+pad_size),pad_size:(y_size+pad_size)] = ori_array.astype(numpy.float64)
    lg_array[0:pad_size,:] = numpy.flipud(lg_array[pad_size:2*pad_size,:])
    lg_array[(x_size+pad_size):(x_size+2*pad_size),:] = numpy.flipud(lg_array[x_size:(x_size+pad_size),:])
    lg_array[:,0:pad_size] = numpy.fliplr(lg_array[:,pad_size:2*pad_size])
    lg_array[:,(y_size+pad_size):(y_size+2*pad_size)] = numpy.fliplr(lg_array[:,y_size:(y_size+pad_size)])
    return lg_array


#
# Classes.
#

def FittingException(Exception):

    def __init__(self, msg):
        Exception.__init__(self, msg)

        
## PeakFinder
#
# Base class for peak finding. This handles identification of peaks in an image.
#
# If you want to modify this with a custom peak finder or an alternative
# way to estimate the background, the recommended approach is to sub-class this
# class and then override backgroundEstimator() and peakFinder().
#
class PeakFinder(object):
    
    # Hard-wired defaults.
    unconverged_dist = 5.0  # Distance between peaks for marking as unconverged (this is multiplied by parameters.sigma)
    new_peak_dist = 1.0     # Minimum allowed distance between new peaks and current peaks.
    
    ## __init__
    #
    # This is called once at the start of analysis to initialize the
    # parameters that will be used for peak fitting.
    #
    # @param parameters A parameters object.
    # @param margin The amount of margin added to the edge of the image.
    #
    def __init__(self, parameters):

        # member variables.
        self.background = None                                              # Current estimate of the image background.
        self.find_max_radius = 5                                            # Radius (in pixels) over which the maxima is maximal.
        self.image = None                                                   # The original image.
        self.iterations = parameters.iterations                             # Maximum number of cycles of peak finding, fitting and subtraction to perform.
        self.margin = PeakFinderFitter.margin                               # Size of the unanalyzed "edge" around the image.
        self.neighborhood = PeakFinder.unconverged_dist * parameters.sigma  # Radius for marking neighbors as unconverged.
        self.new_peak_radius = PeakFinder.new_peak_dist                     # Minimum allowed distance between new peaks and current peaks.
        self.parameters = parameters                                        # Keep access to the parameters object.
        self.peak_mask = None                                               # Mask for limiting peak identification to a particular AOI.
        self.sigma = parameters.sigma                                       # Peak sigma (in pixels).
        self.taken = None                                                   # Spots in the image where a peak has already been added.
        self.threshold = parameters.threshold                               # Peak minimum threshold (height, in camera units).

    ## cleanUp
    #
    def cleanUp(self):
        pass

    ## backgroundEstimator
    #
    # This method does the actual background estimation.
    #
    # Override this if you want to change how the background is estimated.
    #
    def backgroundEstimator(self, image):
        return estimateBackground(image)  # A simple low pass background estimator.
        
    ## findPeaks
    #
    # Finds the peaks in an image & adds to the current list of peaks.
    #
    # @param no_bg_image The current background subtracted image.
    # @param peaks The current list of peaks.
    #
    # @return [True/False if new peaks were added to the current list, the new peaks]
    #
    def findPeaks(self, no_bg_image, peaks):

        # Identify local maxima in the image.
        new_peaks = self.peakFinder(no_bg_image)

        # Fill in initial values for peak height, background and sigma.
        new_peaks = util_c.initializePeaks(new_peaks,         # The new peaks.
                                           self.image,        # The original image.
                                           self.background,   # The current estimate of the background.
                                           self.sigma)        # The starting sigma value.

        # Update new peak identification threshold (if necessary).
        # Also, while threshold is greater than min_threshold we
        # are automatically not done.
        found_new_peaks = False
        if (self.cur_threshold > self.threshold):
            self.cur_threshold -= self.threshold
            found_new_peaks = True

        # If we did not find any new peaks then we may be done.
        if (new_peaks.shape[0] == 0):
            return [found_new_peaks, peaks]

        # Add new peaks to the current list of peaks if it exists,
        # otherwise these peaks become the current list.
        if isinstance(peaks, numpy.ndarray):
            merged_peaks = util_c.mergeNewPeaks(peaks,
                                                new_peaks,
                                                self.new_peak_radius,
                                                self.neighborhood)
        
            # If none of the new peaks are valid then we may be done.
            if (merged_peaks.shape[0] == peaks.shape[0]):
                return [found_new_peaks, merged_peaks]
            else:
                return [True, merged_peaks]
        else:
            return [True, new_peaks]

    ## newImage
    #
    # This is called once at the start of the analysis of a new image.
    #
    # @param image A 2D numpy array.
    #
    def newImage(self, new_image):

        # Make a copy of the starting image.
        self.image = numpy.copy(new_image)
        
        # Reset taken mask.
        self.taken = numpy.zeros(new_image.shape, dtype=numpy.int32) 

        # Initialize new peak minimum threshold.
        if(self.iterations>4):
            self.cur_threshold = 4.0 * self.threshold
        else:
            self.cur_threshold = float(self.iterations) * self.threshold

        # Create mask to limit peak finding to a user defined sub-region of the image.
        if self.peak_mask is None:
            parameters = self.parameters
            self.peak_mask = numpy.ones(new_image.shape)
            if hasattr(parameters, "x_start"):
                self.peak_mask[0:parameters.x_start,:] = 0.0
            if hasattr(parameters, "x_stop"):
                self.peak_mask[parameters.x_stop:-1,:] = 0.0
            if hasattr(parameters, "y_start"):
                self.peak_mask[:,0:parameters.y_start] = 0.0
            if hasattr(parameters, "y_stop"):
                self.peak_mask[:,parameters.y_stop:-1] = 0.0

    ## peakFinder
    #
    # This method does the actual peak finding.
    #
    # Override this if you want to change the peak finding behaviour.
    #
    def peakFinder(self, no_bg_image):

        # Mask the image so that peaks are only found in the AOI.
        masked_image = no_bg_image * self.peak_mask
        
        # Identify local maxima in the masked image.
        [new_peaks, self.taken] = util_c.findLocalMaxima(masked_image,
                                                         self.taken,
                                                         self.cur_threshold,
                                                         self.find_max_radius,
                                                         self.margin)
        return new_peaks

    ## subtractBackground
    #
    # This subtracts a smoothed background out of an image. As a side
    # effect it also updates self.background.
    #
    # @param image The image to subtract the background from.
    #
    # @returns The image with the background subtracted.
    #
    def subtractBackground(self, image):
        self.background = self.backgroundEstimator(image)
        return image - self.background


## PeakFitter
#
# Base class for peak fitting. This handles refinement of the parameters of
# the peaks that were identified with PeakFinder.
#
# Override the peakFitter function if you want to
# change how peak finding is performed.
#
class PeakFitter(object):

    ## __init__
    #
    # @param parameters A (fitting) parameters object.
    #
    def __init__(self, parameters):

        self.image = None        # The image for peak fitting.
        self.scmos_cal = None    # sCMOS calibration data.

        self.neighborhood = parameters.sigma*PeakFinder.unconverged_dist  # Radius for marking neighbors as unconverged.
        self.sigma = parameters.sigma                                     # Peak sigma (in pixels).
        self.threshold = parameters.threshold                             # Peak minimum threshold (height, in camera units).

        # Initialize Z fitting parameters if necessary.
        if (hasattr(parameters, "model")) and (parameters.model == "Z"):
            wx_params = params.getWidthParams(parameters, "x", for_mu_Zfit = True)
            wy_params = params.getWidthParams(parameters, "y", for_mu_Zfit = True)
            [min_z, max_z] = params.getZRange(parameters)

            if(parameters.orientation == "inverted"):
                initZParams(wx_params, wy_params, min_z, max_z)
            else:
                initZParams(wy_params, wx_params, min_z, max_z)

    ## cleanUp
    #
    def cleanUp(self):
        pass

    ## fitPeaks
    #
    # Performs a single iteration of peak fitting.
    #
    # @param peaks A numpy array of peaks to fit.
    #
    # @return [updated peaks, updated residual]
    #
    def fitPeaks(self, peaks):
            
        # Fit to update peak locations.
        [fit_peaks, residual, iterations] = self.peakFitter(peaks)
        fit_peaks = multi_c.getGoodPeaks(fit_peaks,
                                         0.9*self.threshold,
                                         0.5*self.sigma)
            
        # Remove peaks that are too close to each other & refit.
        fit_peaks = util_c.removeClosePeaks(fit_peaks, self.sigma, self.neighborhood)
        [fit_peaks, residual, iterations] = self.peakFitter(fit_peaks)

        fit_peaks = multi_c.getGoodPeaks(fit_peaks,
                                         0.9 * self.threshold,
                                         0.5 * self.sigma)

        return [fit_peaks, residual]

    ## newImage
    #
    # @param new_image A new image (2D numpy array).
    #
    def newImage(self, new_image):
        self.image = numpy.copy(new_image)

    ## peakFitter
    #
    # This method does the actual peak fitting. It is overridden
    # in the sub-class to do the peak fitting.
    #
    # See for example:
    #   3d_daostorm/find_peaks.py
    #
    def peakFitter(self, peaks):
        pass


## PeakFinderFitter
#
# Base class to encapsulate peak finding and fitting. For this to self.peak_finder
# must be set to a PeakFinder object, and self.peak_fitter must be set to a
# PeakFitter object.
#
# To get an idea of how all the pieces are supposed to go together, please see:
#
#   3d_daostorm/find_peaks.py
#   sCMOS/find_peaks.py
#
class PeakFinderFitter():

    margin = 10   # Size of the unanalyzed edge around the image. This is also
                  #  a constant in the C libraries, so if you change this you
                  #  also need to change that.

    ## __init__
    #
    # @param parameters A parameters object
    #
    def __init__(self, parameters):
        self.iterations = parameters.iterations
        self.peak_finder = False           # A sub-class of PeakFinder.
        self.peak_fitter = False           # A sbu-class of PeakFitter.

    ## analyzeImage
    #
    # @param image The image to analyze.
    # @param save_residual (Optional) Save the residual image after peak fitting, default is False.
    #
    # @return [Found peaks, Image residual]
    #
    def analyzeImage(self, new_image, save_residual = False, verbose = False):
        [image, residual] = self.newImage(new_image)

        self.peak_finder.newImage(image)
        self.peak_fitter.newImage(image)

        if save_residual:
            resid_dax = daxwriter.DaxWriter("residual.dax",
                                            residual.shape[0],
                                            residual.shape[1])

        peaks = False
        for i in range(self.iterations):
            if save_residual:
                resid_dax.addFrame(residual)

            no_bg_image = self.peak_finder.subtractBackground(residual)
            [found_new_peaks, peaks] = self.peak_finder.findPeaks(no_bg_image, peaks)
            if isinstance(peaks, numpy.ndarray):
                [peaks, residual] = self.peak_fitter.fitPeaks(peaks)

            if verbose:
                if isinstance(peaks, numpy.ndarray):
                    print " peaks:", i, found_new_peaks, peaks.shape[0]
                else:
                    print " peaks:", i, found_new_peaks, "NA"

            if not found_new_peaks:
                break

        if save_residual:
            resid_dax.addFrame(residual)
            resid_dax.close()

        if isinstance(peaks, numpy.ndarray):
            peaks[:,util_c.getXCenterIndex()] -= float(self.margin)
            peaks[:,util_c.getYCenterIndex()] -= float(self.margin)

        return [peaks, residual]

    ## cleanUp
    #
    # Clean up at the end of analysis.
    #
    def cleanUp(self):
        self.peak_finder.cleanUp()
        self.peak_fitter.cleanUp()

    ## getConvergedPeaks
    #
    # @param peaks A 1D numpy array containing the peaks.
    #
    # @return A 1D numpy array containing only the converged peaks.
    #
    def getConvergedPeaks(self, peaks, verbose = False):
        if (peaks.shape[0] > 0):
            status_index = util_c.getStatusIndex()
            mask = (peaks[:,status_index] == 1.0)  # 0.0 = running, 1.0 = converged.
            if verbose:
                print " ", numpy.sum(mask), "converged out of", peaks.shape[0]
            return peaks[mask,:]
        else:
            return peaks

    ## newImage
    #
    def newImage(self, new_image):
        lg_image = padArray(new_image, self.margin)
        return [lg_image, lg_image]


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
