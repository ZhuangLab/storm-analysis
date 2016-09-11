#!/usr/bin/python
#
# 02/11
#
# Simple Python interface to multi_fit.c for
# ease of performance testing.
#
# 07/11
#
# Modified to reflect changes in multi_fit.c
#
# 09/16
#
# Remove static variables from the C library.
#
# Hazen
#

import ctype
import math
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.loadclib as loadclib

multi = loadclib.loadCLibrary(os.path.dirname(__file__), "multi_fit")


# C interface definition
multi.getResidual.argtypes = [ctypes.c_void_p,
                              ndpointer(dtype=numpy.float64)]
multi.getResults.argtypes = [ctypes.c_void_p,
                             ndpointer(dtype=numpy.float64)]
multi.getUnconverged.restype = [ctypes.c_void_p,
                                ctypes.c_int]
multi.initialize.argtypes = [ndpointer(dtype=numpy.float64),
                             ndpointer(dtype=numpy.float64),
                             ctypes.c_int,
                             ctypes.c_int]
multi.initialize.restype = ctypes.c_void_p

multi.initializeZParameters.argtypes = [ctypes.c_void_p,
                                        ndpointer(dtype=numpy.float64), 
                                        ndpointer(dtype=numpy.float64),
                                        ctypes.c_double,
                                        ctypes.c_double]
multi.iterate2DFixed.argtypes = [ctypes.c_void_p]
multi.iterate2D.argtypes = [ctypes.c_void_p]
multi.iterate3D.argtypes = [ctypes.c_void_p]
multi.iterateZ.argtypes = [ctypes.c_void_p]
multi.newImage.argtypes = [ctypes.c_void_p,
                           ndpointer(dtype=numpy.float64)]


# Globals
default_tol = 1.0e-6

peakpar_size = utilC.getNPeakPar()


class MultiFitterException(Exception):
    
    def __init__(self, message):
        Exception.__init__(self, message)


class MultiFitter(object):
    """
    This is designed to be used as follows:

    (1) At the start of the analysis, create a single instance of the appropriate fitting sub-class.
    (2) For each new image, call newImage() once.
    (3) For each iteration of peak fittings, call doFit() with peaks that you want to fit to the image.
    (4) After calling doFit() you can remove peaks that did not fit will with getGoodPeaks().
    (5) After calling doFit() you can use getResidual() to get the current image minus the fit peaks.
    (6) Call cleanup() when you are done with this object and plan to throw it away.

    As all the static variables have been removed from the C library you should 
    be able to use several of these objects simultaneuosly for fitting.
    """
    def __init__(self, scmos_cal = None, wx_params = None, wy_params = None, min_z = None, max_z = None, verbose = False)

        self.max_z = max_z
        self.min_z = min_z
        self.scmos_cal = scmos_cal
        self.verbose = verbose
        self.wx_params = wx_params
        self.wy_params = wy_params

        self.mfit = None

        # Default clamp parameters.
        #
        # These set the scale for how much these parameters
        # can change in a single fitting iteration.
        #
        self.clamp = numpy.array([1000.0,   # Height
                                  1.0,      # x position
                                  0.3,      # width in x
                                  1.0,      # y position
                                  0.3,      # width in y
                                  100.0,    # background
                                  0.1])     # z position

    def cleanup(self):
        multi.cleanup(self.mfit)
        self.mfit = None

    def doFit(self, peaks, max_iterations = 200):

        # Initialize C library with new peaks.
        multi.newPeaks(self.mfit,
                       peaks,
                       numpy.ascontiguousarray(peaks.shape[0]))

        # Iterate fittings.
        i = 0
        self.iterate()
        while(multi.getUnconverged(self.mfit) and (i < max_iterations)):
            if self.verbose and ((i%20)==0):
                print("iteration", i)
            self.iterate()
            i += 1

        if self.verbose:
            if (i == max_iterations):
                print(" Failed to converge in:", i, multi.getUnconverged(self.mfit))
            else:
                print(" Multi-fit converged in:", i, multi.getUnconverged(self.mfit))
            print("")

        # Get updated peak values back from the C library.
        fit_peaks = numpy.ascontiguousarray(numpy.zeros(peaks.shape))
        multi.getResults(self.mfit, fit_peaks)
        return fit_peaks

    def getGoodPeaks(self, peaks, min_height, min_width):
        """
        Create a new list from peaks containing only those peaks that meet 
        the specified criteria for minimum peak height and width.
        """
        if(peaks.shape[0]>0):
            min_width = 0.5 * min_width

            status_index = utilC.getStatusIndex()
            height_index = utilC.getHeightIndex()
            xwidth_index = utilC.getXWidthIndex()
            ywidth_index = utilC.getYWidthIndex()

            if verbose:
                tmp = numpy.ones(peaks.shape[0])                
                print("getGoodPeaks")
                for i in range(peaks.shape[0]):
                    print(i, peaks[i,0], peaks[i,1], peaks[i,3], peaks[i,2], peaks[i,4])
                print("Total peaks:", numpy.sum(tmp))
                print("  fit error:", numpy.sum(tmp[(peaks[:,status_index] != 2.0)]))
                print("  min height:", numpy.sum(tmp[(peaks[:,height_index] > min_height)]))
                print("  min width:", numpy.sum(tmp[(peaks[:,xwidth_index] > min_width) & (peaks[:,ywidth_index] > min_width)]))
                print("")
            mask = (peaks[:,7] != status_index) & (peaks[:,height_index] > min_height) & (peaks[:,xwidth_index] > min_width) & (peaks[:,ywidth_index] > min_width)
            return peaks[mask,:]
        else:
            return peaks

    def getResidual(self):
        residual = numpy.ascontiguousarray(numpy.zeros(self.scmos_cal.shape))
        multi.getResidual(self.mfit, residual)
        return residual
    
    def iterate(self):
        """
        Sub classes override this to use the correct fitting function.
        """
        pass

    def initMFit(self, image);
        """
        This initializes the C fitting library. You can call this directly, but
        the idea is that it will get called automatically the first time that you
        provide a new image for fitting.
        """.
        if self.scmos_cal is None:
            self.scmos_cal = numpy.ascontiguousarray(numpy.zeros(image.shape))0202
        else:
            self.scmos_cal = numpy.ascontiguousarray(self.scmos_cal)

        if (image.shape[0] != self.scmos_cal.shape[0]) or (image.shape[1] != self.scmos_cal.shape[1]):
            raise MultiFitterException("Image shape and sCMOS calibration shape do not match.")

        self.mfit = multi.initialize(self.scmos_cal,
                                     numpy.ascontiguousarray(self.clamp),
                                     self.scmos_cal.shape[1],
                                     self.scmos_cal.shape[0])

        if self.wx_params is not None:
            multi.initializeZParameters(self.mfit,
                                        numpy.ascontiguousarray(self.wx_params),
                                        numpy.ascontiguousarray(self.wy_params),
                                        self.min_z,
                                        self.max_z)

    def newImage(self, image):

        if self.mfit is None:
            self.initMFit(image)

        multi.newImage(self.mfit, image)


class MultiFitter2DFixed(MultiFitter):
    """
    Fit with a fixed peak width.
    """
    def iterate(self):
        multi.iterate2DFixed(self.mfit)


class MultiFitter2D(MultiFitter):
    """
    Fit with a variable peak width (of the same size in X and Y).
    """
    def iterate(self):
        multi.iterate2D(self.mfit)

        
class MultiFitter3D(MultiFitter):
    """
    Fit with peak width that can change independently in X and Y.
    """
    def iterate(self):
        multi.iterate3D(self.mfit)

        
class MultiFitterZ(MultiFitter):
    """
    Fit with peak width that varies in X and Y as a function of Z.
    """
    def iterate(self):
        multi.iterateZ(self.mfit)


#
# Other functions.
#
def calcSxSy(wx_params, wy_params, z):
    """
    Return sigma x and sigma y given the z calibration parameters.
    """
    zx = (z - wx_params[1])/wx_params[2]
    sx = 0.5 * wx_params[0] * math.sqrt(1.0 + zx*zx + wx_params[3]*zx*zx*zx + wx_params[4]*zx*zx*zx*zx)
    zy = (z - wy_params[1])/wy_params[2]
    sy = 0.5 * wy_params[0] * math.sqrt(1.0 + zy*zy + wy_params[3]*zy*zy*zy + wy_params[4]*zy*zy*zy*zy)
    return [sx, sy]



### TO REMOVE ####

def fitStats(results):
    total = results.shape[0]
    bad = numpy.count_nonzero(results[:,n_params] == 2.0)
    good = numpy.count_nonzero(results[:,n_params] == 1.0)
    unc = numpy.count_nonzero(results[:,n_params] == 0.0)
    return [good, bad, unc, total]

def getConvergedPeaks(peaks, min_height, min_width):
    if(peaks.shape[0]>0):
        min_width = 0.5 * min_width
        mask = (peaks[:,7] == 1.0) & (peaks[:,0] > min_height) & (peaks[:,2] > min_width) & (peaks[:,4] > min_width)
        return peaks[mask,:]
    else:
        return peaks

def getGoodPeaks(peaks, min_height, min_width, verbose = False):
    if(peaks.shape[0]>0):
        min_width = 0.5 * min_width
        if verbose:
            tmp = numpy.ones(peaks.shape[0])
            print("getGoodPeaks")
            for i in range(peaks.shape[0]):
                print(i, peaks[i,0], peaks[i,1], peaks[i,3], peaks[i,2], peaks[i,4])
            print("Total peaks:", numpy.sum(tmp))
            print("  fit error:", numpy.sum(tmp[(peaks[:,7] != 2.0)]))
            print("  min height:", numpy.sum(tmp[(peaks[:,0] > min_height)]))
            print("  min width:", numpy.sum(tmp[(peaks[:,2] > min_width) & (peaks[:,4] > min_width)]))
            print("")
        mask = (peaks[:,7] != 2.0) & (peaks[:,0] > min_height) & (peaks[:,2] > min_width) & (peaks[:,4] > min_width)
        return peaks[mask,:]
    else:
        return peaks

def getResidual(im_size_x, im_size_y):
    arr = numpy.ascontiguousarray(numpy.zeros(im_size_x*im_size_y))
    multi.getResidual(arr)
    return numpy.reshape(arr, (im_size_x, im_size_y))

def getResults(size):
    arr1 = numpy.ascontiguousarray(numpy.zeros(size*resultspar_size))
    multi.getResults(arr1)
    return numpy.reshape(arr1, (-1, resultspar_size))

# Pretty prints peak fitting results.
def prettyPrint(results):
    peaks = results
    n_res = peaks.shape[0]
    names = ["H", "XC", "SX", "YC", "SY", "BG", "Z", "ST", "ERR"]
    for i in range(n_res):
        f = "good"
        if(peaks[i,7] == 0.0):
            f = "unconverged"
        if(peaks[i,7] == 2.0):
            f = "error"
        print("Peak", i, f)
        for j, name in enumerate(names):
            print("   ", name, "%.3f" % (peaks[i,j]))


# Print fitting stats.
def printStats(results):
    [good, bad, unc, total] = fitStats(results)
    print("%d g: %.3f b: %.3f u: %.3f" % (total, float(good)/total, float(bad)/total, float(unc)/total))


##
# Fitting Functions.
##

# Generic fitting function.
#
def _doFit_(fitfn, data, scmos_cal, peaks, tolerance, max_iters, verbose, zfit):
    if verbose:
        print("_doFit_")
        for i in range(peaks.shape[0]):
            print(i, peaks[i,0], peaks[i,1], peaks[i,5])
        print("")

    if isinstance(scmos_cal, numpy.ndarray):
        c_scmos_cal = numpy.ascontiguousarray(scmos_cal)
    else:
        c_scmos_cal = numpy.ascontiguousarray(numpy.zeros(data.shape))
    n_peaks = int(peaks.size/resultspar_size)
    c_data = numpy.ascontiguousarray(data)
    c_peaks = numpy.ascontiguousarray(peaks)
    if verbose:
        print("initializing, ", n_peaks, "peaks")
    multi.initialize(c_data,
                     c_scmos_cal,
                     c_peaks,
                     tolerance,
                     c_data.shape[1],
                     c_data.shape[0], 
                     n_peaks,
                     zfit)

    if verbose:
        print("fitting")
    i = 1
    fitfn()
    while(multi.getUnconverged() and (i < max_iters)):
        if verbose and ((i%20)==0):
            print("iteration", i)
        #if (i<10):
        #    print "iteration", i
        #    print "  ", i, multi.getUnconverged()
        fitfn()
        i += 1

    #if verbose:
    if verbose:
        if (i==max_iters):
            print(" Failed to converge in:", i, multi.getUnconverged())
        else:
            print(" Multi-fit converged in:", i, multi.getUnconverged())
        print("")

    fit = getResults(n_peaks)
    res = getResidual(data.shape[0], data.shape[1])
    multi.cleanup()
    return fit, res, i


# Fits a single gaussian peak of fixed width.
# 
# Assumed centered in the image as this is mostly designed to be
# used for testing purposes.
def fitSingleGaussian2DFixed(data, sigma, scmos_cal = False, tolerance = default_tol, verbose = False):
    peaks = [numpy.max(data)-numpy.min(data),
             0.5 * data.shape[0],
             sigma,
             0.5 * data.shape[1],
             sigma,
             numpy.min(data),
             0,
             0,
             0.0]
    return _doFit_(multi.iterate2DFixed, data, scmos_cal, numpy.array(peaks), tolerance, 100, verbose, 0)

# Fits a single gaussian peak with the same width in x and y.
# 
# Assumed centered in the image as this is mostly designed to be
# used for testing purposes.
def fitSingleGaussian2D(data, sigma, scmos_cal = False, tolerance = default_tol, verbose = False):
    peaks = [numpy.max(data)-numpy.min(data),
             0.5 * data.shape[0],
             sigma,
             0.5 * data.shape[1],
             sigma,
             numpy.min(data),
             0,
             0,
             0.0]
    return _doFit_(multi.iterate2D, data, scmos_cal, numpy.array(peaks), tolerance, 100, verbose, 0)

# Fits a single gaussian peak w/ varying width in x and y.
# 
# Assumed centered in the image as this is mostly designed to be
# used for testing purposes.
def fitSingleGaussian3D(data, sigma, scmos_cal = False, tolerance = default_tol, verbose = False):
    peaks = [numpy.max(data)-numpy.min(data),
             0.5 * data.shape[0],
             sigma,
             0.5 * data.shape[1],
             sigma,
             numpy.min(data),
             0,
             0,
             0.0]
    return _doFit_(multi.iterate3D, data, scmos_cal, numpy.array(peaks), tolerance, 100, verbose, 0)

# Fits a single gaussian peak w/ x, y width depending on z.
# 
# Assumed centered in the image as this is mostly designed to be
# used for testing purposes.
def fitSingleGaussianZ(data, wx_params, wy_params, scmos_cal = False, tolerance = default_tol, verbose = False):
    c_wx = numpy.ascontiguousarray(wx_params)
    c_wy = numpy.ascontiguousarray(wy_params)
    multi.initializeZParameters(wx_params, wy_params)

    [sx, sy] = calcSxSy(wx_params, wy_params, 0.0)
    peaks = [numpy.max(data)-numpy.min(data),
             0.5 * data.shape[0],
             sx,
             0.5 * data.shape[1],
             sy,
             numpy.min(data),
             0.0,
             0,
             0.0]
    return _doFit_(multi.iterateZ, data, scmos_cal, numpy.array(peaks), tolerance, 100, verbose, 1)

# Fits multiple gaussian peaks of fixed width.
#
def fitMultiGaussian2DFixed(data, peaks, scmos_cal, tolerance = default_tol, max_iters = 200, verbose = False):
    return _doFit_(multi.iterate2DFixed, data, scmos_cal, peaks, tolerance, max_iters, verbose, 0)

# Fits multiple gaussian peaks w/ varying width
# but symmetric in x and y.
#
def fitMultiGaussian2D(data, peaks, scmos_cal, tolerance = default_tol, max_iters = 200, verbose = False):
    return _doFit_(multi.iterate2D, data, scmos_cal, peaks, tolerance, max_iters, verbose, 0)

# Fits multiple gaussian peaks w/ varying width in x and y.
#
def fitMultiGaussian3D(data, peaks, scmos_cal, tolerance = default_tol, max_iters = 200, verbose = False):
    return _doFit_(multi.iterate3D, data, scmos_cal, peaks, tolerance, max_iters, verbose, 0)

# Fits multiple gaussian peaks w/ x, y width varying based on z.
#
def fitMultiGaussianZ(data, peaks, scmos_cal, tolerance = default_tol, max_iters = 200, verbose = False):
    return _doFit_(multi.iterateZ, data, scmos_cal, peaks, tolerance, max_iters, verbose, 1)


###
# Z fit setup functions.
###

# Initialize parameters for Z fitting.
def initZParams(wx_params, wy_params, min_z, max_z):
    if(min_z > max_z):
        print("Bad z range detected", min_z, max_z, "switching to defaults.")
        min_z = -0.5
        max_z = 0.5
    c_wx = numpy.ascontiguousarray(wx_params)
    c_wy = numpy.ascontiguousarray(wy_params)
    multi.initializeZParameters(c_wx, c_wy, min_z, max_z)


#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
