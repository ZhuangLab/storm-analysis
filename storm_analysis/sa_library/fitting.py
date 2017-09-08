#!/usr/bin/env python
"""
This contains some of what is common to all of the peak finding algorithms.

Hazen 03/17
"""

import numpy
import os
import scipy
import scipy.ndimage

import storm_analysis.sa_library.daxwriter as daxwriter
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.parameters as params


#
# Functions.
#
def estimateBackground(image, size = 8):
    """
    A simple background estimator.

    image - A 2D numpy array containing the image.
    size - (Optional) The size of the filter to use, default is 8.

    returns - A 2D numpy array containing the image background estimate.
    """
    background = scipy.ndimage.filters.gaussian_filter(image, (size, size))
    return background

def loadSCMOSData(calibration_filename, margin):
    """
    Load camera calibration data.
    
    Note: Gain is expected to be in units of ADU per photo-electron.
    """
    [offset, variance, gain] = numpy.load(calibration_filename)

    # Pad out camera calibration data to the final image size.
    lg_offset = padArray(offset, margin)
    lg_variance = padArray(variance, margin)
    lg_gain = padArray(gain, margin)

    return [lg_offset, lg_variance, lg_gain]

def padArray(ori_array, pad_size):
    """
    Pads out an array to a large size.

    ori_array - A 2D numpy array.
    pad_size - The number of elements to add to each of the "sides" of the array.
    
    The padded 2D numpy array.
    """
    if (pad_size > 0):
        [x_size, y_size] = ori_array.shape
        lg_array = numpy.ones((x_size+2*pad_size,y_size+2*pad_size))
        lg_array[pad_size:(x_size+pad_size),pad_size:(y_size+pad_size)] = ori_array.astype(numpy.float64)
        lg_array[0:pad_size,:] = numpy.flipud(lg_array[pad_size:2*pad_size,:])
        lg_array[(x_size+pad_size):(x_size+2*pad_size),:] = numpy.flipud(lg_array[x_size:(x_size+pad_size),:])
        lg_array[:,0:pad_size] = numpy.fliplr(lg_array[:,pad_size:2*pad_size])
        lg_array[:,(y_size+pad_size):(y_size+2*pad_size)] = numpy.fliplr(lg_array[:,y_size:(y_size+pad_size)])
        return lg_array
    
    else:
        return ori_array


#
# Classes.
#
def FittingException(Exception):

    def __init__(self, msg):
        Exception.__init__(self, msg)

        
class PeakFinder(object):
    """
    Base class for peak finding. This handles identification of peaks in an image.

    If you want to modify this with a custom peak finder or an alternative
    way to estimate the background, the recommended approach is to sub-class this
    class and then override backgroundEstimator() and peakFinder().
    """
    # Hard-wired defaults.
    
    unconverged_dist = 5.0  # Distance between peaks for marking as unconverged (this is multiplied by parameters.sigma)
    new_peak_dist = 1.0     # Minimum allowed distance between new peaks and current peaks.
    
    def __init__(self, parameters):
        """
        This is called once at the start of analysis to initialize the
        parameters that will be used for peak fitting.
 
        parameters - A parameters object.
        """
        # Initialized from parameters.
        self.find_max_radius = parameters.getAttr("find_max_radius", 5)     # Radius (in pixels) over which the maxima is maximal.
        self.iterations = parameters.getAttr("iterations")                  # Maximum number of cycles of peak finding, fitting and subtraction to perform.
        self.sigma = parameters.getAttr("sigma")                            # Peak sigma (in pixels).
        self.threshold = parameters.getAttr("threshold")                    # Peak minimum threshold (height, in camera units).
        self.z_value = parameters.getAttr("z_value", 0.0)                   # The starting z value to use for peak fitting.

        # Other member variables.
        self.background = None                                              # Current estimate of the image background.
        self.image = None                                                   # The original image.
        self.margin = PeakFinderFitter.margin                               # Size of the unanalyzed "edge" around the image.
        self.neighborhood = PeakFinder.unconverged_dist * self.sigma        # Radius for marking neighbors as unconverged.
        self.new_peak_radius = PeakFinder.new_peak_dist                     # Minimum allowed distance between new peaks and current peaks.
        self.parameters = parameters                                        # Keep access to the parameters object.
        self.peak_locations = None                                          # Initial peak locations, as explained below.
        self.peak_mask = None                                               # Mask for limiting peak identification to a particular AOI.
        self.taken = None                                                   # Spots in the image where a peak has already been added.

        
        #
        # This is for is you already know where your want fitting to happen, as
        # for example in a bead calibration movie and you just want to use the
        # approximate locations as inputs for fitting.
        #
        # peak_locations is a text file with the peak x, y, height and background
        # values as white spaced columns (x and y positions are in pixels as
        # determined using visualizer).
        #
        # 1.0 2.0 1000.0 100.0
        # 10.0 5.0 2000.0 200.0
        # ...
        #
        if parameters.hasAttr("peak_locations"):

            peak_filename = parameters.getAttr("peak_locations")
            if os.path.exists(peak_filename):
                print("Using peak starting locations specified in", peak_filename)
            elif os.path.exists(os.path.basename(peak_filename)):
                peak_filename = os.path.basename(peak_filename)
                print("Using peak starting locations specified in", peak_filename)

            # Only do one cycle of peak finding as we'll always return the same locations.
            if (self.iterations != 1):
                print("WARNING: setting number of iterations to 1!")
                self.iterations = 1

            # Load peak x,y locations.
            peak_locs = numpy.loadtxt(peak_filename, ndmin = 2)
            print("Loaded", peak_locs.shape[0], "peak locations")

            # Create peak array.
            self.peak_locations = numpy.zeros((peak_locs.shape[0],
                                               utilC.getNPeakPar()))
            self.peak_locations[:,utilC.getXCenterIndex()] = peak_locs[:,1] + self.margin
            self.peak_locations[:,utilC.getYCenterIndex()] = peak_locs[:,0] + self.margin
            self.peak_locations[:,utilC.getHeightIndex()] = peak_locs[:,2]
            self.peak_locations[:,utilC.getBackgroundIndex()] = peak_locs[:,3]

            self.peak_locations[:,utilC.getXWidthIndex()] = numpy.ones(peak_locs.shape[0]) * self.sigma
            self.peak_locations[:,utilC.getYWidthIndex()] = numpy.ones(peak_locs.shape[0]) * self.sigma

    def backgroundEstimator(self, image):
        """
        This method does the actual background estimation.

        Override this if you want to change how the background is estimated.
        """
        return estimateBackground(image)  # A simple low pass background estimator.

    def cleanUp(self):
        pass

    def findPeaks(self, no_bg_image, peaks):
        """
        Finds the peaks in an image & adds to the current list of peaks.
   
        no_bg_image - The current background subtracted image.
        peaks - The current list of peaks.
    
        return - [True/False if new peaks were added to the current list, the new peaks]
        """

        # Use pre-specified peak locations if available, e.g. bead calibration.
        if self.peak_locations is not None:
            new_peaks = self.peak_locations
            
        # Otherwise, identify local maxima in the image and initialize fitting parameters.
        else:
            new_peaks = self.peakFinder(no_bg_image)

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
            merged_peaks = self.mergeNewPeaks(peaks, new_peaks)
        
            # If none of the new peaks are valid then we may be done.
            if (merged_peaks.shape[0] == peaks.shape[0]):
                return [found_new_peaks, merged_peaks]
            else:
                return [True, merged_peaks]
        else:
            return [True, new_peaks]

    def mergeNewPeaks(self, peaks, new_peaks):
        """
        Merge new peaks into the current list of peaks.
        """
        return utilC.mergeNewPeaks(peaks,
                                   new_peaks,
                                   self.new_peak_radius,
                                   self.neighborhood)
        
    def newImage(self, new_image):
        """
        This is called once at the start of the analysis of a new image.
        
        new_image - A 2D numpy array.
        """

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
            self.peak_mask = numpy.ones(new_image.shape)
            if self.parameters.hasAttr("x_start"):
                self.peak_mask[0:self.parameters.getAttr("x_start")+self.margin,:] = 0.0
            if self.parameters.hasAttr("x_stop"):
                self.peak_mask[self.parameters.getAttr("x_stop")+self.margin:-1,:] = 0.0
            if self.parameters.hasAttr("y_start"):
                self.peak_mask[:,0:self.parameters.getAttr("y_start")+self.margin] = 0.0
            if self.parameters.hasAttr("y_stop"):
                self.peak_mask[:,self.parameters.getAttr("y_stop")+self.margin:-1] = 0.0

    def peakFinder(self, no_bg_image):
        """
        This method does the actual peak finding.
        
        Override this if you want to change the peak finding behaviour.
        """
        # Mask the image so that peaks are only found in the AOI.
        masked_image = no_bg_image * self.peak_mask
        
        # Identify local maxima in the masked image.
        [new_peaks, self.taken] = utilC.findLocalMaxima(masked_image,
                                                        self.taken,
                                                        self.cur_threshold,
                                                        self.find_max_radius,
                                                        self.margin)

        # Fill in initial values for peak height, background and sigma.
        new_peaks = utilC.initializePeaks(new_peaks,         # The new peaks.
                                          self.image,        # The original image.
                                          self.background,   # The current estimate of the background.
                                          self.sigma,        # The starting sigma value.
                                          self.z_value)      # The starting z value.
            
        return new_peaks

    def subtractBackground(self, image, bg_estimate):
        """
        This subtracts a smoothed background out of an image. As a side
        effect it also updates self.background.
        
        image - The image to subtract the background from.
        bg_estimate - An estimate of the background.

        return - The image with the background subtracted.
        """

        # If we are provided with an estimate of the background
        # then just use it.
        if bg_estimate is not None:
            self.background = bg_estimate

        # Otherwise make our own estimate.
        else:
            self.background = self.backgroundEstimator(image)
            
        return image - self.background


class PeakFitter(object):
    """
    Base class for peak fitting. This handles refinement of the
    parameters of the peaks that were identified with PeakFinder.

    The actual fitting is done by an the self.mfitter object, this
    is primarily just a wrapper for the self.mfitter object.
    """

    def __init__(self, parameters):
        """
        parameters - A (fitting) parameters object.
        """

        self.image = None        # The image for peak fitting.
        self.mfitter = None      # An instance of a sub-class of the MultiFitter class.
        self.scmos_cal = None    # sCMOS calibration data.

        # Z fitting parameters.
        self.max_z = None
        self.min_z = None
        self.wx_params = None
        self.wy_params = None

        self.sigma = parameters.getAttr("sigma")                          # Peak sigma (in pixels).
        self.threshold = parameters.getAttr("threshold")                  # Peak minimum threshold (height, in camera units).

        self.neighborhood = self.sigma*PeakFinder.unconverged_dist        # Radius for marking neighbors as unconverged.

        # Initialize Z fitting parameters if necessary.
        if (parameters.getAttr("model", "na") == "Z"):
            [self.wx_params, self.wy_params] = parameters.getWidthParams(for_mu_Zfit = True)
            [self.min_z, self.max_z] = parameters.getZRange()

    def cleanUp(self):
        self.mfitter.cleanup()

    def fitPeaks(self, peaks):
        """
        Performs a single iteration of peak fitting.
        
        peaks - A numpy array of peaks to fit.
    
        return - [updated peaks, updated residual]
        """
        
        # Fit to update peak locations.
        [fit_peaks, residual] = self.peakFitter(peaks)
        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks,
                                              0.9*self.threshold,
                                              0.5*self.sigma)
        
        # Remove peaks that are too close to each other & refit.
        fit_peaks = utilC.removeClosePeaks(fit_peaks, self.sigma, self.neighborhood)
        [fit_peaks, residual] = self.peakFitter(fit_peaks)

        fit_peaks = self.mfitter.getGoodPeaks(fit_peaks,
                                              0.9 * self.threshold,
                                              0.5 * self.sigma)
        
        return [fit_peaks, residual]

    def newImage(self, new_image):
        """
        new_image - A new image (2D numpy array).
        """
        self.mfitter.newImage(new_image)

    def peakFitter(self, peaks):
        """
        This method does the actual peak fitting.
        """
        fit_peaks = self.mfitter.doFit(peaks)
        residual = self.mfitter.getResidual()
        return [fit_peaks, residual]


class PeakFinderFitter():
    """
    Base class to encapsulate peak finding and fitting. For this to self.peak_finder
    must be set to a PeakFinder object, and self.peak_fitter must be set to a
    PeakFitter object.

    To get an idea of how all the pieces are supposed to go together, please see:
      3d_daostorm/find_peaks.py
      sCMOS/find_peaks.py
    """
    
    margin = 10   # Size of the unanalyzed edge around the image. This is also
                  #  a constant in the C libraries, so if you change this you
                  #  also need to change that.

    def __init__(self, parameters):
        """
        parameters - A parameters object.
        """
        self.peak_finder = False           # A sub-class of PeakFinder.
        self.peak_fitter = False           # A sub-class of PeakFitter.

    def analyzeImage(self, movie_reader, save_residual = False, verbose = False):
        """
        movie_reader - std_analysis.MovieReader object.
        save_residual - (Optional) Save the residual image after peak fitting, default is False.

        return - [Found peaks, Image residual]
        """
        bg_estimate = movie_reader.getBackground()
        new_image = movie_reader.getFrame()
        
        #
        # Pad out arrays so that we can better analyze localizations
        # near the edge of the original image.
        #
        image = padArray(new_image, self.margin)
        residual = padArray(new_image, self.margin)
        if bg_estimate is not None:
            bg_estimate = padArray(bg_estimate, self.margin)

        self.peak_finder.newImage(image)
        self.peak_fitter.newImage(image)

        if save_residual:
            resid_dax = daxwriter.DaxWriter("residual.dax",
                                            residual.shape[0],
                                            residual.shape[1])

        peaks = False
        for i in range(self.peak_finder.iterations):
            if save_residual:
                resid_dax.addFrame(residual)

            no_bg_image = self.peak_finder.subtractBackground(residual, bg_estimate)
            [found_new_peaks, peaks] = self.peak_finder.findPeaks(no_bg_image, peaks)
            if isinstance(peaks, numpy.ndarray):
                [peaks, residual] = self.peak_fitter.fitPeaks(peaks)

            if verbose:
                if isinstance(peaks, numpy.ndarray):
                    print(" peaks:", i, found_new_peaks, peaks.shape[0])
                else:
                    print(" peaks:", i, found_new_peaks, "NA")

            if not found_new_peaks:
                break

        if save_residual:
            resid_dax.addFrame(residual)
            resid_dax.close()

        if isinstance(peaks, numpy.ndarray):
            peaks[:,utilC.getXCenterIndex()] -= float(self.margin)
            peaks[:,utilC.getYCenterIndex()] -= float(self.margin)

        return [peaks, residual]

    def cleanUp(self):
        self.peak_finder.cleanUp()
        self.peak_fitter.cleanUp()

    def getConvergedPeaks(self, peaks, verbose = False):
        """
        peaks - A 1D numpy array containing the peaks.
        
        return - A 1D numpy array containing only the converged peaks.
        """
        if (peaks.shape[0] > 0):
            status_index = utilC.getStatusIndex()
            mask = (peaks[:,status_index] == 1.0)  # 0.0 = running, 1.0 = converged.
            if verbose:
                print(" ", numpy.sum(mask), "converged out of", peaks.shape[0])
            return peaks[mask,:]
        else:
            return peaks


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
