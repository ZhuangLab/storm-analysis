#!/usr/bin/env python
"""
Python interface to the dao_fit C library.

Hazen
"""

import ctypes
import math
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.loadclib as loadclib


#
# The Python definitions of the C structures in sa_library/multi_fit.h
#
class fitData(ctypes.Structure):
    _fields_ = [('n_dposv', ctypes.c_int),
                ('n_iterations', ctypes.c_int),
                ('n_margin', ctypes.c_int),
                ('n_neg_fi', ctypes.c_int),
                ('n_neg_height', ctypes.c_int),
                ('n_neg_width', ctypes.c_int),

                ('margin', ctypes.c_int),
                ('nfit', ctypes.c_int),
                ('image_size_x', ctypes.c_int),
                ('image_size_y', ctypes.c_int),

                ('xoff', ctypes.c_double),
                ('yoff', ctypes.c_double),
                ('zoff', ctypes.c_double),
                
                ('tolerance', ctypes.c_double),
                
                ('bg_counts', ctypes.POINTER(ctypes.c_int)),
                
                ('bg_data', ctypes.POINTER(ctypes.c_double)),
                ('f_data', ctypes.POINTER(ctypes.c_double)),
                ('scmos_term', ctypes.POINTER(ctypes.c_double)),
                ('x_data', ctypes.POINTER(ctypes.c_double)),

                ('clamp_start', (ctypes.c_double*7)),
                
                ('fit', ctypes.c_void_p),
                ('fit_model', ctypes.c_void_p)]


def loadDaoFitC():
    daofit = loadclib.loadCLibrary("storm_analysis.sa_library", "dao_fit")
    
    # These are from sa_library/multi_fit.c
    daofit.mFitGetFitImage.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]
    
    daofit.mFitGetResidual.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]
    
    daofit.mFitGetResults.argtypes = [ctypes.c_void_p,
                                      ndpointer(dtype=numpy.float64)]

    daofit.mFitGetUnconverged.argtypes = [ctypes.c_void_p]
    daofit.mFitGetUnconverged.restype = ctypes.c_int

    daofit.mFitNewImage.argtypes = [ctypes.c_void_p,
                                    ndpointer(dtype=numpy.float64)]

    # These are from sa_library/dao_fit.c
    daofit.cleanup.argtypes = [ctypes.c_void_p]
        
    daofit.initialize.argtypes = [ndpointer(dtype=numpy.float64),
                                  ndpointer(dtype=numpy.float64),
                                  ctypes.c_double,
                                  ctypes.c_int,
                                  ctypes.c_int]
    daofit.initialize.restype = ctypes.POINTER(fitData)
    #daofit.initialize.restype = ctypes.c_void_p
    
    daofit.initializeZParameters.argtypes = [ctypes.c_void_p,
                                             ndpointer(dtype=numpy.float64), 
                                             ndpointer(dtype=numpy.float64),
                                             ctypes.c_double,
                                             ctypes.c_double]

    daofit.iterate2DFixed.argtypes = [ctypes.POINTER(fitData)]
    
    daofit.iterate2D.argtypes = [ctypes.c_void_p]
    
    daofit.iterate3D.argtypes = [ctypes.c_void_p]
    
    daofit.iterateZ.argtypes = [ctypes.c_void_p]
    
    daofit.newPeaks.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype=numpy.float64),
                                ctypes.c_int]

    return daofit
    

def printFittingInfo(mfit, spacing = "  "):
    """
    Print out some of the information the C fitting library keeps track of.
    """
    print(spacing, mfit.contents.n_dposv, "fits lost to Cholesky failure")
    print(spacing, mfit.contents.n_margin, "fits lost to image margin")
    print(spacing, mfit.contents.n_neg_fi, "fits lost to negative value in fit function")
    print(spacing, mfit.contents.n_neg_height, "fits lost to negative height")
    print(spacing, mfit.contents.n_neg_width, "fits lost to negative width")

class MultiFitterException(Exception):
    
    def __init__(self, message):
        Exception.__init__(self, message)


class MultiFitterBase(object):
    """
    Base class to make it easier to share some functionality with Spliner.
    """
    def __init__(self, scmos_cal = None, verbose = False, min_z = None, max_z = None, **kwds):
        super(MultiFitterBase, self).__init__(**kwds)
        self.clib = None
        self.default_tol = 1.0e-6
        self.im_shape = None
        self.iterations = 0
        self.max_z = max_z
        self.mfit = None
        self.min_z = min_z
        self.scmos_cal = scmos_cal
        self.verbose = verbose
        
        # Default clamp parameters.
        #
        # These set the (initial) scale for how much these parameters
        # can change in a single fitting iteration.
        #
        self.clamp = numpy.array([1.0,  # Height (Note: This is relative to the initial guess).
                                  1.0,  # x position
                                  0.3,  # width in x
                                  1.0,  # y position
                                  0.3,  # width in y
                                  1.0,  # background (Note: This is relative to the initial guess).
                                  0.1]) # z position
        
    def cleanup(self, spacing = "  ", verbose = True):
        """
        This just prints the analysis statistics, it does not do any actual cleanup.
        """
        if self.mfit is not None:
            if verbose:
                printFittingInfo(self.mfit, spacing = spacing)
                print(spacing, self.iterations, "fitting iterations.")

    def doFit(self, peaks, max_iterations = 200):
        """
        This is where the fitting actually happens.
        """
        # Initialize C library with new peaks.
        self.newPeaks(peaks)

        # Iterate fittings.
        i = 0
        self.iterate()
        while(self.getUnconverged() and (i < max_iterations)):
            if self.verbose and ((i%20)==0):
                print("iteration", i)
            self.iterate()
            i += 1

        if self.verbose:
            if (i == max_iterations):
                print(" Failed to converge in:", i, self.getUnconverged())
            else:
                print(" Multi-fit converged in:", i, self.getUnconverged())
            print("")

        # Get number of fitting iterations.
        self.getIterations()
        
        # Get updated peak values back from the C library.
        return self.getResults(peaks.shape)

    def getFitImage(self):
        """
        Get the fit image, i.e. f(x), an image created from drawing all of
        the current fits into a 2D array.
        """
        fit_image = numpy.zeros(self.im_shape)
        self.clib.mFitGetFitImage(self.mfit, fit_image)
        return fit_image

    def getIterations(self):
        """
        Update iterations and reset C library counter. The idea anyway
        is that the Python counter won't overflow where as the C counter
        might, particularly on a 32 bit system.
        """
        self.iterations += self.mfit.contents.n_iterations
        self.mfit.contents.n_iterations = 0
        return self.iterations
        
    def getResidual(self):
        """
        Get the residual, the data minus the fit image, xi - f(x).
        """
        residual = numpy.ascontiguousarray(numpy.zeros(self.im_shape))
        self.clib.mFitGetResidual(self.mfit, residual)
        return residual

    def getResults(self, input_peaks_shape):
        """
        Get the peaks, presumably after a few rounds of fitting to improve
        their parameters.
        """
        fit_peaks = numpy.ascontiguousarray(numpy.zeros(input_peaks_shape))
        self.clib.mFitGetResults(self.mfit, fit_peaks)
        return fit_peaks

    def getUnconverged(self):
        """
        Return the number of fits that have not yet converged.
        """
        return self.clib.mFitGetUnconverged(self.mfit)

    def initializeC(self, image):
        """
        This initializes the C fitting library. You can call this directly, but
        the idea is that it will get called automatically the first time that you
        provide a new image for fitting.
        """
        if self.scmos_cal is None:
            self.scmos_cal = numpy.ascontiguousarray(numpy.zeros(image.shape))
        else:
            self.scmos_cal = numpy.ascontiguousarray(self.scmos_cal)

        if (image.shape[0] != self.scmos_cal.shape[0]) or (image.shape[1] != self.scmos_cal.shape[1]):
            raise MultiFitterException("Image shape and sCMOS calibration shape do not match.")

        self.im_shape = self.scmos_cal.shape
        
    def iterate(self):
        """
        Sub classes override this to use the correct fitting function.
        """
        raise MultiFitterException("iterate() method not defined.")

    def newImage(self, image):
        """
        Initialize the C fitter with new data, an image as a numpy.ndarray. If the 
        fitter does not exist thie will also initialize the C fitter.
        """
        if self.mfit is None:
            self.initializeC(image)
        else:
            if (image.shape[0] != self.im_shape[0]) or (image.shape[1] != self.im_shape[1]):
                raise MultiFitterException("Current image shape and the original image shape are not the same.")

        self.clib.mFitNewImage(self.mfit, image)

    def newPeaks(self, peaks):
        """
        Sub classes override this to provide analysis specific peak initialization.
        """
        raise MultiFitterException("newPeaks() method not defined.")
    
    
class MultiFitter(MultiFitterBase):
    """
    This is designed to be used as follows:

    (1) At the start of the analysis, create a single instance of the appropriate fitting sub-class.
    (2) For each new image, call newImage() once.
    (3) For each iteration of peak fittings, call doFit() with peaks that you want to fit to the image.
    (4) After calling doFit() you can remove peaks that did not fit well with getGoodPeaks().
    (5) After calling doFit() you can use getResidual() to get the current image minus the fit peaks.
    (6) Call cleanup() when you are done with this object and plan to throw it away.

    As all the static variables have been removed from the C library you should 
    be able to use several of these objects simultaneuosly for fitting.

    All of the parameters are optional, use None if they are not relevant.
    """
    def __init__(self, wx_params = None, wy_params = None, **kwds):
        super(MultiFitter, self).__init__(**kwds)

        self.wx_params = wx_params
        self.wy_params = wy_params

        self.clib = loadDaoFitC()

    def cleanup(self, verbose = True):
        super(MultiFitter, self).cleanup(verbose = verbose)
        if self.mfit is not None:
            self.clib.cleanup(self.mfit)
            self.mfit = None

    def getGoodPeaks(self, peaks,  min_width):
        """
        Create a new list from peaks containing only those peaks that meet 
        the specified criteria for minimum peak width.

        FIXME: Using sigma threshold we don't really have a value for minimum
               height. Do we need to restore a height filter?
        """
        if(peaks.shape[0]>0):
            min_width = 0.5 * min_width

            status_index = utilC.getStatusIndex()
            height_index = utilC.getHeightIndex()
            xwidth_index = utilC.getXWidthIndex()
            ywidth_index = utilC.getYWidthIndex()

            if self.verbose:
                tmp = numpy.ones(peaks.shape[0])                
                print("getGoodPeaks")
                for i in range(peaks.shape[0]):
                    print(i, peaks[i,0], peaks[i,1], peaks[i,3], peaks[i,2], peaks[i,4], peaks[i,7])
                print("Total peaks:", numpy.sum(tmp))
                print("  fit error:", numpy.sum(tmp[(peaks[:,status_index] != 2.0)]))
                print("  min width:", numpy.sum(tmp[(peaks[:,xwidth_index] > min_width) & (peaks[:,ywidth_index] > min_width)]))
                print("")
            mask = (peaks[:,status_index] != 2.0) & (peaks[:,xwidth_index] > min_width) & (peaks[:,ywidth_index] > min_width)
            return peaks[mask,:]
        else:
            return peaks

    def initializeC(self, image):
        """
        This initializes the C fitting library. You can call this directly, but
        the idea is that it will get called automatically the first time that you
        provide a new image for fitting.
        """
        super(MultiFitter, self).initializeC(image)
        
        self.mfit = self.clib.initialize(self.scmos_cal,
                                         numpy.ascontiguousarray(self.clamp),
                                         self.default_tol,
                                         self.scmos_cal.shape[1],
                                         self.scmos_cal.shape[0])

        if self.wx_params is not None:
            self.clib.initializeZParameters(self.mfit,
                                            numpy.ascontiguousarray(self.wx_params),
                                            numpy.ascontiguousarray(self.wy_params),
                                            self.min_z,
                                            self.max_z)

    def newPeaks(self, peaks):
        """
        Pass new peaks to the C library.
        """
        self.clib.newPeaks(self.mfit,
                           numpy.ascontiguousarray(peaks),
                           peaks.shape[0])


class MultiFitter2DFixed(MultiFitter):
    """
    Fit with a fixed peak width.
    """
    def iterate(self):
        self.clib.iterate2DFixed(self.mfit)


class MultiFitter2D(MultiFitter):
    """
    Fit with a variable peak width (of the same size in X and Y).
    """
    def iterate(self):
        self.clib.iterate2D(self.mfit)

        
class MultiFitter3D(MultiFitter):
    """
    Fit with peak width that can change independently in X and Y.
    """
    def iterate(self):
        self.clib.iterate3D(self.mfit)

        
class MultiFitterZ(MultiFitter):
    """
    Fit with peak width that varies in X and Y as a function of Z.
    """
    def iterate(self):
        self.clib.iterateZ(self.mfit)


#
# Other functions.
#
# FIXME: This is cruft.
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
