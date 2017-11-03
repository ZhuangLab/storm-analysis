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
                ('n_non_decr', ctypes.c_int),

                ('jac_size', ctypes.c_int),
                ('margin', ctypes.c_int),
                ('max_nfit', ctypes.c_int),
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

                ('working_peak', ctypes.c_void_p),
                ('fit', ctypes.c_void_p),

                ('fit_model', ctypes.c_void_p),

                ('fn_add_peak', ctypes.c_void_p),
                ('fn_alloc_peaks', ctypes.c_void_p),
                ('fn_calc_JH', ctypes.c_void_p),
                ('fn_calc_peak_shape', ctypes.c_void_p),
                ('fn_check', ctypes.c_void_p),
                ('fn_copy_peak', ctypes.c_void_p),
                ('fn_free_peaks', ctypes.c_void_p),                
                ('fn_subtract_peak', ctypes.c_void_p),
                ('fn_update', ctypes.c_void_p)]


def loadDaoFitC():
    daofit = loadclib.loadCLibrary("storm_analysis.sa_library", "dao_fit")
    
    # These are from sa_library/multi_fit.c
    daofit.mFitGetFitImage.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]
    
    daofit.mFitGetNError.argtypes = [ctypes.c_void_p]
    daofit.mFitGetNError.restype = ctypes.c_int
    
    daofit.mFitGetPeakPropertyDouble.argtypes = [ctypes.c_void_p,
                                                 ndpointer(dtype=numpy.float64),
                                                 ctypes.c_char_p]

    daofit.mFitGetPeakPropertyInt.argtypes = [ctypes.c_void_p,
                                              ndpointer(dtype=numpy.int32),
                                              ctypes.c_char_p]
    
    daofit.mFitGetResidual.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]
    
    daofit.mFitGetUnconverged.argtypes = [ctypes.c_void_p]
    daofit.mFitGetUnconverged.restype = ctypes.c_int

    daofit.mFitIterateLM.argtypes = [ctypes.c_void_p]
    daofit.mFitIterateOriginal.argtypes = [ctypes.c_void_p]

    daofit.mFitNewBackground.argtypes = [ctypes.c_void_p,
                                         ndpointer(dtype=numpy.float64)]
    
    daofit.mFitNewImage.argtypes = [ctypes.c_void_p,
                                    ndpointer(dtype=numpy.float64)]

    daofit.mFitRemoveErrorPeaks.argtypes = [ctypes.c_void_p]
                                             
    daofit.mFitSetPeakStatus.argtypes = [ctypes.c_void_p,
                                         ndpointer(dtype=numpy.int32)]

    # These are from sa_library/dao_fit.c
    daofit.daoCleanup.argtypes = [ctypes.c_void_p]
        
    daofit.daoInitialize.argtypes = [ndpointer(dtype=numpy.float64),
                                     ndpointer(dtype=numpy.float64),
                                     ctypes.c_double,
                                     ctypes.c_int,
                                     ctypes.c_int]
    daofit.daoInitialize.restype = ctypes.POINTER(fitData)

    daofit.daoInitialize2DFixed.argtypes = [ctypes.c_void_p]
    daofit.daoInitialize2D.argtypes = [ctypes.c_void_p]
    daofit.daoInitialize3D.argtypes = [ctypes.c_void_p]
    daofit.daoInitializeZ.argtypes = [ctypes.c_void_p]
    
    daofit.daoInitializeZ.argtypes = [ctypes.c_void_p,
                                      ndpointer(dtype=numpy.float64), 
                                      ndpointer(dtype=numpy.float64),
                                      ctypes.c_double,
                                      ctypes.c_double]
    
    daofit.daoNewPeaks.argtypes = [ctypes.c_void_p,
                                   ndpointer(dtype=numpy.float64),
                                   ctypes.c_char_p,
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
    print(spacing, mfit.contents.n_non_decr, "fits reset due to non-decreasing error (LM).")


class MultiFitterException(Exception):
    pass


class MultiFitterBase(object):
    """
    Base class for fitting multiple possibly overlapping localizations.
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
        self.peak_properties = {"background" : "float",
                                "error" : "float",
                                "height" : "float",
                                "status" : "int",
                                "x" : "float",
                                "xwidth" : "float",
                                "y" : "float",
                                "ywidth" : "float",
                                "z" : "float"}
        
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

    def doFit(self, max_iterations = 200):
        """
        This is where the fitting actually happens.

        FIXME: Why we always do at least one iteration? I guess it doesn't
               matter because this should be a NOP as the C library will
               not do anything if all the peaks have converged.
        """
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

    def getFitImage(self):
        """
        Get the fit image, i.e. f(x), an image created from drawing all of
        the current fits into a 2D array.
        """
        fit_image = numpy.ascontiguousarray(numpy.zeros(self.im_shape, dtype = numpy.float64))
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

    def getNError(self):
        """
        Return the number of peaks in the C library that are in the ERROR state.
        """
        return self.clib.mFitGetNError(self.mfit)
    
    def getNFit(self):
        """
        Return the current number of peaks that the C library is handling.
        """
        return self.mfit.contents.nfit

    def getNFitMax(self):
        """
        Return the current maximum number of peaks. 
        
        Note this is not a fixed value as the C library can dynamically
        increase this. This method is primarily for testing purposes.
        """
        return self.mfit.contents.max_nfit

    def getPeakProperty(self, p_name):
        """
        Return a numpy array containing the requested property.
        """
        if not p_name in self.properties:
            raise MultiFitterException("No such property '" + p_name + "'")

        if(self.properties[p_name] == "float"):
            values = numpy.ascontiguousarray(numpy.zeros(self.getNFit(), dtype = numpy.float64))
            return self.clib.mFitGetPeakPropertyDouble(self.mfit,
                                                       values,
                                                       ctypes.c_char_p(p_name.encode()))
        elif(self.properties[p_name] == "int"):
            values = numpy.ascontiguousarray(numpy.zeros(self.getNFit(), dtype = numpy.int32))
            return self.clib.mFitGetPeakPropertyInt(self.mfit,
                                                    values,
                                                    ctypes.c_char_p(p_name.encode()))

    def getResidual(self):
        """
        Get the residual, the data minus the fit image, xi - f(x).
        """
        residual = numpy.ascontiguousarray(numpy.zeros(self.im_shape, dtype = numpy.float64))
        self.clib.mFitGetResidual(self.mfit, residual)
        return residual

#    def getResults(self, input_peaks_shape):
#        """
#        Get the peaks, presumably after a few rounds of fitting to improve
#        their parameters.
#        """
#        fit_peaks = numpy.ascontiguousarray(numpy.zeros(input_peaks_shape))
#        self.clib.mFitGetResults(self.mfit, fit_peaks)
#        return fit_peaks

    def getUnconverged(self):
        """
        Return the number of fits that have not yet converged.
        """
        return self.clib.mFitGetUnconverged(self.mfit)

    def initializeC(self, image):
        """
        This initializes the C fitting library.

        It needs the image in order to know what size arrays to create
        as we won't always have SCMOS calibration data.
        """
        if self.scmos_cal is None:
            self.scmos_cal = numpy.ascontiguousarray(numpy.zeros(image.shape), dtype = numpy.float64)
        else:
            self.scmos_cal = numpy.ascontiguousarray(self.scmos_cal, dtype = numpy.float64)

        if (image.shape[0] != self.scmos_cal.shape[0]) or (image.shape[1] != self.scmos_cal.shape[1]):
            raise MultiFitterException("Image shape and sCMOS calibration shape do not match.")

        self.im_shape = self.scmos_cal.shape

    def isInitialized(self):
        return (self.mfit != None)
    
    def iterate(self):
        """
        Sub classes override this to use the correct fitting function.
        """
        raise MultiFitterException("iterate() method not defined.")

    def newBackground(self, background):
        """
        Update the current background estimate.
        """
        if (background.shape[0] != self.im_shape[0]) or (background.shape[1] != self.im_shape[1]):
            raise MultiFitterException("Background image shape and the original image shape are not the same.")
        self.clib.mFitNewBackground(self.mfit,
                                    numpy.ascontiguousarray(background, dtype = numpy.float64))
        
    def newImage(self, image):
        """
        Initialize C fitter with a new image.
        """
        if (image.shape[0] != self.im_shape[0]) or (image.shape[1] != self.im_shape[1]):
            raise MultiFitterException("Current image shape and the original image shape are not the same.")

        self.clib.mFitNewImage(self.mfit,
                               numpy.ascontiguousarray(image, dtype = numpy.float64))

    def newPeaks(self, peaks, peaks_type):
        """
        Sub classes override this to provide analysis specific peak addition.
        """
        raise MultiFitterException("newPeaks() method not defined.")

    def removeErrorPeaks(self):
        """
        Instruct the C library to remove all the peaks in the ERROR state
        from the list of peaks that it is maintaining.
        """
        self.clib.mFitRemoveErrorPeaks(self.mfit)
        
    def setPeakStatus(self, status):
        """
        Set the status (RUNNING, CONVERGED, ERROR) of the peaks in the C library.
        """
        assert (status.size == self.getNFit())
        self.clib.mFitSetPeakStatus(self.mfit,
                                    numpy.ascontiguousarray(status, dtype = numpy.int32))
            

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
            self.clib.daoCleanup(self.mfit)
            self.mfit = None

    def initializeC(self, image):
        """
        This initializes the C fitting library.
        """
        super(MultiFitter, self).initializeC(image)
        
        self.mfit = self.clib.daoInitialize(self.scmos_cal,
                                            numpy.ascontiguousarray(self.clamp, dtype = numpy.float64),
                                            self.default_tol,
                                            self.scmos_cal.shape[1],
                                            self.scmos_cal.shape[0])

    def iterate(self):
        self.clib.mFitIterateLM(self.mfit)
        #self.clib.mFitIterateOriginal(self.mfit)

    def newPeaks(self, peaks, peaks_type):
        """
        Pass new peaks to add to the C library.
        """
        self.clib.daoNewPeaks(self.mfit,
                              numpy.ascontiguousarray(peaks),
                              ctypes.c_char_p(peaks_type.encode()),
                              peaks.shape[0])

    def rescaleZ(self, peaks):
        """
        Convert Z from fitting units to microns.
        """
        return peaks


class MultiFitter2DFixed(MultiFitter):
    """
    Fit with a fixed peak width.
    """
    def initializeC(self, image):
        super(MultiFitter2DFixed, self).initializeC(image)
        self.clib.daoInitialize2DFixed(self.mfit)


class MultiFitter2D(MultiFitter):
    """
    Fit with a variable peak width (of the same size in X and Y).
    """
    def initializeC(self, image):
        super(MultiFitter2D, self).initializeC(image)
        self.clib.daoInitialize2D(self.mfit)
        
        
class MultiFitter3D(MultiFitter):
    """
    Fit with peak width that can change independently in X and Y.
    """
    def initializeC(self, image):
        super(MultiFitter3D, self).initializeC(image)
        self.clib.daoInitialize3D(self.mfit)
        
        
class MultiFitterZ(MultiFitter):
    """
    Fit with peak width that varies in X and Y as a function of Z.
    """
    def initializeC(self, image):
        super(MultiFitterZ, self).initializeC(image)
        self.clib.daoInitializeZ(self.mfit,
                                 numpy.ascontiguousarray(self.wx_params),
                                 numpy.ascontiguousarray(self.wy_params),
                                 self.min_z,
                                 self.max_z)


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
