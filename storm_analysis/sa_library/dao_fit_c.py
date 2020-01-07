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
                ('n_lost', ctypes.c_int),
                ('n_margin', ctypes.c_int),
                ('n_neg_fi', ctypes.c_int),
                ('n_neg_height', ctypes.c_int),
                ('n_non_converged', ctypes.c_int),
                ('n_non_decr', ctypes.c_int),

                ('fit_size_x', ctypes.c_int),
                ('fit_size_y', ctypes.c_int),
                ('image_size_x', ctypes.c_int),
                ('image_size_y', ctypes.c_int),                
                ('jac_size', ctypes.c_int),
                ('max_nfit', ctypes.c_int),
                ('nfit', ctypes.c_int),
                ('roi_n_index', ctypes.c_int),

                ('min_height', ctypes.c_double),
                
                ('xoff', ctypes.c_double),
                ('yoff', ctypes.c_double),
                ('zoff', ctypes.c_double),
                
                ('tolerance', ctypes.c_double),
                
                ('bg_counts', ctypes.POINTER(ctypes.c_int)),
                ('roi_x_index', ctypes.POINTER(ctypes.c_int)),
                ('roi_y_index', ctypes.POINTER(ctypes.c_int)),
                ('stale', ctypes.POINTER(ctypes.c_int)),

                ('as_xi', ctypes.POINTER(ctypes.c_double)),
                ('bg_data', ctypes.POINTER(ctypes.c_double)),
                ('bg_estimate', ctypes.POINTER(ctypes.c_double)),
                ('err_i', ctypes.POINTER(ctypes.c_double)),
                ('f_data', ctypes.POINTER(ctypes.c_double)),
                ('rqe', ctypes.POINTER(ctypes.c_double)),
                ('scmos_term', ctypes.POINTER(ctypes.c_double)),
                ('t_fi', ctypes.POINTER(ctypes.c_double)),
                ('x_data', ctypes.POINTER(ctypes.c_double)),

                ('working_peak', ctypes.c_void_p),
                ('fit', ctypes.c_void_p),

                ('fit_model', ctypes.c_void_p),

                ('fn_alloc_peaks', ctypes.c_void_p),
                ('fn_calc_JH', ctypes.c_void_p),
                ('fn_calc_peak_shape', ctypes.c_void_p),
                ('fn_check', ctypes.c_void_p),
                ('fn_copy_peak', ctypes.c_void_p),
                ('fn_error_fn', ctypes.c_void_p),
                ('fn_free_peaks', ctypes.c_void_p),
                ('fn_peak_sum', ctypes.c_void_p),
                ('fn_update', ctypes.c_void_p)]


def formatPeaksArbitraryPSF(peaks, peaks_type):
    """
    Input peaks array formatter for arbitrary PSFs.

    Based on peaks_type, create a properly formatted ndarray to pass
    to the C library. This is primarily for internal use by newPeaks().
    """
    # These come from the finder, or the unit test code, create peaks
    # as (N,3) with columns x, y, z.
    #
    # Note: "testing" is designed specifically for use by the unit
    #       tests. If you use this then initial height estimation
    #       for the peaks is not performed.
    #
    if (peaks_type == "testing") or (peaks_type == "finder"):
        c_peaks = numpy.stack((peaks["x"],
                               peaks["y"],
                               peaks["z"]), axis = 1)

    # These come from pre-specified peak fitting locations, create peaks
    # as (N,5) with columns x, y, z, background, height.
    #
    elif (peaks_type == "text") or (peaks_type == "hdf5"):
        c_peaks = numpy.stack((peaks["x"],
                               peaks["y"],
                               peaks["z"],
                               peaks["background"],
                               peaks["height"]), axis = 1)
    else:
        raise MultiFitterException("Unknown peaks type '" + peaks_type + "'")

    return numpy.ascontiguousarray(c_peaks, dtype = numpy.float64)


def formatPeaksGaussianPSF(peaks, peaks_type):
    """
    Input peaks array formatter for Gaussian PSFs.

    Based on peaks_type, create a properly formatted ndarray to pass
    to the C library. This is primarily for internal use by newPeaks().
    """
    # These come from the finder, or the unit test code, create peaks
    # as (N,4) with columns x, y, z, and sigma.
    #
    if (peaks_type == "testing") or (peaks_type == "finder"):
        c_peaks = numpy.stack((peaks["x"],
                               peaks["y"],
                               peaks["z"],
                               peaks["sigma"]), axis = 1)
        
    # These come from pre-specified peak fitting locations, create peaks
    # as (N,7) with columns x, y, z, background, height, xsigma, ysigma.
    #
    elif (peaks_type == "text") or (peaks_type == "hdf5"):
        c_peaks = numpy.stack((peaks["x"],
                               peaks["y"],
                               peaks["z"],
                               peaks["background"],
                               peaks["height"],
                               peaks["xsigma"],
                               peaks["ysigma"]), axis = 1)
    else:
        raise MultiFitterException("Unknown peaks type '" + peaks_type + "'")

    return numpy.ascontiguousarray(c_peaks, dtype = numpy.float64)
    

def loadDaoFitC():
    daofit = loadclib.loadCLibrary("dao_fit")
    
    # These are from sa_library/multi_fit.c
    daofit.mFitAnscombeTransformImage.argtypes = [ctypes.c_void_p]

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

    daofit.mFitNewBackground.argtypes = [ctypes.c_void_p,
                                         ndpointer(dtype=numpy.float64)]
    
    daofit.mFitNewImage.argtypes = [ctypes.c_void_p,
                                    ndpointer(dtype=numpy.float64)]

    daofit.mFitRemoveErrorPeaks.argtypes = [ctypes.c_void_p]

    daofit.mFitRemoveRunningPeaks.argtypes = [ctypes.c_void_p]

    daofit.mFitSetPeakStatus.argtypes = [ctypes.c_void_p,
                                         ndpointer(dtype=numpy.int32)]

    # These are from sa_library/dao_fit.c
    daofit.daoCleanup.argtypes = [ctypes.c_void_p]
        
    daofit.daoInitialize.argtypes = [ndpointer(dtype=numpy.float64),
                                     ndpointer(dtype=numpy.float64),
                                     ctypes.c_double,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_int]
    daofit.daoInitialize.restype = ctypes.POINTER(fitData)

    daofit.daoInitialize2DFixed.argtypes = [ctypes.c_void_p]    
    daofit.daoInitialize2DFixedALS.argtypes = [ctypes.c_void_p]
    daofit.daoInitialize2DFixedLS.argtypes = [ctypes.c_void_p]
    daofit.daoInitialize2DFixedDWLS.argtypes = [ctypes.c_void_p]
    daofit.daoInitialize2DFixedFWLS.argtypes = [ctypes.c_void_p]
    
    daofit.daoInitialize2D.argtypes = [ctypes.c_void_p,
                                       ctypes.c_double,
                                       ctypes.c_double]
    daofit.daoInitialize2DALS.argtypes = [ctypes.c_void_p,
                                          ctypes.c_double,
                                          ctypes.c_double]
    daofit.daoInitialize2DLS.argtypes = [ctypes.c_void_p,
                                         ctypes.c_double,
                                         ctypes.c_double]
    daofit.daoInitialize2DDWLS.argtypes = [ctypes.c_void_p,
                                           ctypes.c_double,
                                           ctypes.c_double]
    daofit.daoInitialize2DFWLS.argtypes = [ctypes.c_void_p,
                                           ctypes.c_double,
                                           ctypes.c_double]
     
    daofit.daoInitialize3D.argtypes = [ctypes.c_void_p,
                                       ctypes.c_double,
                                       ctypes.c_double]
    
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
    print(spacing, mfit.contents.n_dposv, "fits reset due to Cholesky failure.")
    print(spacing, mfit.contents.n_margin, "fits reset due to image margin.")
    print(spacing, mfit.contents.n_neg_fi, "fits reset due to negative value in fit function.")
    print(spacing, mfit.contents.n_neg_height, "fits reset due to negative height.")
    print(spacing, mfit.contents.n_non_decr, "fits reset due to non-decreasing error (LM).")
    print(spacing, mfit.contents.n_non_converged, "fits did not converge.")
    print(spacing, mfit.contents.n_lost, "fits were lost.")

class MultiFitterException(Exception):
    pass


class MultiFitter(object):
    """
    Base class for fitting multiple possibly overlapping localizations. This is designed to be 
    used as follows:

     1. At the start of the analysis, create a single instance of the appropriate fitting sub-class.
     2. For each new image, call newImage() once.
     3. Provide an estimate of the background with newBackground().
     4. Add peaks to fit with newPeaks().
     5. Call doFit() to fit the peaks.
     6. Use multiple calls to getPeakProperties() to get the properties you are interested in.
     7. Call cleanup() when you are done with this object and plan to throw it away.

    See sa_library/fitting.py for a typical work flow.

    As all the static variables have been removed from the C library you should 
    be able to use several of these objects simultaneuosly for fitting.

    All of the parameters are optional, use None if they are not relevant.
    """
    def __init__(self,
                 rqe = None,
                 scmos_cal = None,
                 verbose = False,
                 min_z = None,
                 max_z = None,
                 **kwds):
        super(MultiFitter, self).__init__(**kwds)
        self.clib = None
        self.default_tol = 1.0e-6
        self.im_shape = None
        self.iterations = 0
        self.max_z = max_z
        self.mfit = None
        self.min_z = min_z
        self.n_proximity = 0
        self.n_significance = 0

        # These are all the peak (localization) properties that the C libraries
        # estimate. Not all C libraries will provide estimates for all of these
        # properties. It is used by the getPeakProperty() method to check that
        # the requested property is valid.
        #
        self.peak_properties = {"background" : "float",
                                "bg_sum" : "float",
                                "error" : "float",
                                "fg_sum" : "float",
                                "height" : "float",
                                "iterations" : "int",
                                "jacobian" : "float",
                                "significance" : "compound",
                                "sum" : "float",
                                "status" : "int",
                                "x" : "float",
                                "xsigma" : "float",
                                "y" : "float",
                                "ysigma" : "float",
                                "z" : "float"}

        self.rqe = rqe
        self.scmos_cal = scmos_cal

        self.verbose = verbose

    def cleanup(self, spacing = "  ", verbose = True):
        """
        This just prints the analysis statistics, it does not do any actual cleanup.
        """
        if self.mfit is not None:
            if verbose:
                printFittingInfo(self.mfit, spacing = spacing)
                print(spacing, self.n_proximity, "peaks lost to proximity filter.")
                print(spacing, self.n_significance, "peaks lost to low significance.")
                print(spacing, self.iterations, "fitting iterations.")

    def doFit(self, max_iterations = 2000):
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
        if not p_name in self.peak_properties:
            raise MultiFitterException("No such property '" + p_name + "'")

        # Properties that are calculated from other properties.
        if(self.peak_properties[p_name] == "compound"):

            # Return 0 length array if there are no localizations.
            if(self.getNFit() == 0):
                return numpy.zeros(0, dtype = numpy.float64)
                
            # Peak significance calculation.
            if(p_name == "significance"):
                bg_sum = self.getPeakProperty("bg_sum")
                fg_sum = self.getPeakProperty("fg_sum")
                return fg_sum/numpy.sqrt(bg_sum)
    
        # Floating point properties.
        elif(self.peak_properties[p_name] == "float"):
            if (p_name == "jacobian"):
                values = numpy.ascontiguousarray(numpy.zeros((self.getNFit(), self.mfit.contents.jac_size),
                                                             dtype = numpy.float64))
                self.clib.mFitGetPeakPropertyDouble(self.mfit,
                                                    values,
                                                    ctypes.c_char_p(p_name.encode()))
            else:
                values = numpy.ascontiguousarray(numpy.zeros(self.getNFit(), dtype = numpy.float64))
                self.clib.mFitGetPeakPropertyDouble(self.mfit,
                                                    values,
                                                    ctypes.c_char_p(p_name.encode()))
            return values

        # Integer properties.
        elif(self.peak_properties[p_name] == "int"):
            values = numpy.ascontiguousarray(numpy.zeros(self.getNFit(), dtype = numpy.int32))
            self.clib.mFitGetPeakPropertyInt(self.mfit,
                                             values,
                                             ctypes.c_char_p(p_name.encode()))
            return values
        
    def getResidual(self):
        """
        Get the residual, the data minus the fit image, xi - f(x).
        """
        residual = numpy.ascontiguousarray(numpy.zeros(self.im_shape, dtype = numpy.float64))
        self.clib.mFitGetResidual(self.mfit, residual)
        return residual

    def getUnconverged(self):
        """
        Return the number of fits that have not yet converged.
        """
        return self.clib.mFitGetUnconverged(self.mfit)

    def incProximityCounter(self, n_inc):
        self.n_proximity += n_inc

    def incSignificanceCounter(self, n_inc):
        self.n_significance += n_inc

    def initializeC(self, image):
        """
        This initializes the C fitting library.

        It needs the image in order to know what size arrays to create
        as we won't always have SCMOS calibration data.
        """
        if self.scmos_cal is None:
            if self.verbose:
                print("Using zeros for sCMOS calibration data.")
            self.scmos_cal = numpy.ascontiguousarray(numpy.zeros(image.shape), dtype = numpy.float64)
        else:
            self.scmos_cal = numpy.ascontiguousarray(self.scmos_cal, dtype = numpy.float64)

        if self.rqe is None:
            if self.verbose:
                print("Using ones for relative quantum efficiency data.")
            self.rqe = numpy.ascontiguousarray(numpy.ones(image.shape), dtype = numpy.float64)
        else:
            self.rqe = numpy.ascontiguousarray(self.rqe, dtype = numpy.float64)

        if (image.shape[0] != self.scmos_cal.shape[0]) or (image.shape[1] != self.scmos_cal.shape[1]):
            raise MultiFitterException("Image shape and sCMOS calibration shape do not match.")

        if (image.shape[0] != self.rqe.shape[0]) or (image.shape[1] != self.rqe.shape[1]):
            raise MultiFitterException("Image shape and RQE shape do not match.")
        
        self.im_shape = self.scmos_cal.shape

    def isInitialized(self):
        return (self.mfit != None)
    
    def iterate(self):
        self.clib.mFitIterateLM(self.mfit)

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
        Instruct the C library to discard all the peaks in the ERROR state
        from the list of peaks that it is maintaining.
        """
        self.clib.mFitRemoveErrorPeaks(self.mfit)

    def removeRunningPeaks(self):
        """
        Instruct the C library to discard all the peaks in the RUNNING state
        from the list of peaks that it is maintaining. This is usually called
        at the end of the analysis after all of the peaks in the ERROR state
        have been removed.
        """
        self.clib.mFitRemoveRunningPeaks(self.mfit)

    def rescaleZ(self, z):
        """
        Convert Z from fitting units to microns.
        """
        return z
        
    def setPeakStatus(self, status):
        """
        Set the status (RUNNING, CONVERGED, ERROR) of the peaks in the C library.
        """
        assert (status.size == self.getNFit())
        self.clib.mFitSetPeakStatus(self.mfit,
                                    numpy.ascontiguousarray(status, dtype = numpy.int32))
            

class MultiFitterArbitraryPSF(MultiFitter):
    """
    Base class for arbitrary PSF fitters (Spliner, PupilFn, PSFFFT)
    """
    def formatPeaks(self, peaks, peaks_type):
        return formatPeaksArbitraryPSF(peaks, peaks_type)
    
    
class MultiFitterGaussian(MultiFitter):
    """
    Base class for Gaussian fitters (3D-DAOSTORM and sCMOS).
    """
    def __init__(self, roi_size = 10, wx_params = None, wy_params = None, **kwds):
        super(MultiFitterGaussian, self).__init__(**kwds)

        self.roi_size = roi_size
        self.wx_params = wx_params
        self.wy_params = wy_params

        self.clib = loadDaoFitC()

    def cleanup(self, verbose = True):
        super(MultiFitterGaussian, self).cleanup(verbose = verbose)
        if self.mfit is not None:
            self.clib.daoCleanup(self.mfit)
            self.mfit = None

    def formatPeaks(self, peaks, peaks_type):
        return formatPeaksGaussianPSF(peaks, peaks_type)
            
    def initializeC(self, image):
        """
        This initializes the C fitting library.
        """
        super(MultiFitterGaussian, self).initializeC(image)
        
        self.mfit = self.clib.daoInitialize(self.rqe,
                                            self.scmos_cal,
                                            self.default_tol,
                                            self.scmos_cal.shape[1],
                                            self.scmos_cal.shape[0],
                                            self.roi_size)

    def newPeaks(self, peaks, peaks_type):
        """
        Pass new peaks to add to the C library.
        """
        c_peaks = self.formatPeaks(peaks, peaks_type)
        self.clib.daoNewPeaks(self.mfit,
                              c_peaks,
                              ctypes.c_char_p(peaks_type.encode()),
                              c_peaks.shape[0])


class MultiFitter2DFixed(MultiFitterGaussian):
    """
    Fit with a fixed peak width.
    """
    def initializeC(self, image):
        super(MultiFitter2DFixed, self).initializeC(image)
        self.clib.daoInitialize2DFixed(self.mfit)

        
class MultiFitter2DFixedALS(MultiFitterGaussian):
    """
    Fit with a fixed peak width using the Anscombe least squares fitting
    error model.
    """
    def initializeC(self, image):
        super(MultiFitter2DFixedALS, self).initializeC(image)
        self.clib.daoInitialize2DFixedALS(self.mfit)

    def newImage(self, image):
         super(MultiFitter2DFixedALS, self).newImage(image)
         self.clib.mFitAnscombeTransformImage(self.mfit)

  
class MultiFitter2DFixedLS(MultiFitterGaussian):
    """
    Fit with a fixed peak width using the Least squares fitting
    error model.
    """
    def initializeC(self, image):
        super(MultiFitter2DFixedLS, self).initializeC(image)
        self.clib.daoInitialize2DFixedLS(self.mfit)


class MultiFitter2DFixedDWLS(MultiFitterGaussian):
    """
    Fit with a fixed peak width using the data weighted least squares 
    fitting error model.
    """
    def initializeC(self, image):
        super(MultiFitter2DFixedDWLS, self).initializeC(image)
        self.clib.daoInitialize2DFixedDWLS(self.mfit)


class MultiFitter2DFixedFWLS(MultiFitterGaussian):
    """
    Fit with a fixed peak width using the fit weighted least squares 
    fitting error model.
    """
    def initializeC(self, image):
        super(MultiFitter2DFixedFWLS, self).initializeC(image)
        self.clib.daoInitialize2DFixedFWLS(self.mfit)

        
class MultiFitter2DFixedNC(MultiFitter2DFixed):
    """
    Fit with a fixed peak width, but without correcting for RQE. More
    specifically we set the RQE correction to 1.0 so that the fitter
    will use the same RQE correction approach as the finder (the
    original image is divided by the RQE). At least in theory this
    will be slightly worse as the statistics will no longer be
    exactly Poisson. In practice it appears that the differences are
    somewhere in the 4th or 5th digit, so pretty small.

    This is primarily for testing.
    """
    def __init__(self, **kwds):
        super(MultiFitter2DFixedNC, self).__init__(**kwds)
        self.rqe = None


class MultiFitter2DFixedDWLSNC(MultiFitter2DFixedDWLS):
    """
    Fit with a fixed peak width using the data weighted least squares 
    fitting error model, but without correcting for RQE. More
    specifically we set the RQE correction to 1.0 so that the fitter
    will use the same RQE correction approach as the finder (the
    original image is divided by the RQE).

    Using this we can test the performance of the combination of mean
    gain and flat field correction with weighted least squares fitting
    as reported by other labs. We'll use the same value for all of the
    gains and use the RQE term as the flat field correction.

    This is primarily for testing.
    """
    def __init__(self, **kwds):
        super(MultiFitter2DFixedDWLSNC, self).__init__(**kwds)
        self.rqe = None
        
        
class MultiFitter2D(MultiFitterGaussian):
    """
    Fit with a variable peak width (of the same size in X and Y).
    """
    def __init__(self, sigma_range = None, **kwds):
        super(MultiFitter2D, self).__init__(**kwds)
        
        self.sigma_range = sigma_range
            
    def initializeC(self, image):
        super(MultiFitter2D, self).initializeC(image)

        width_max = 1.0/(2.0 * self.sigma_range[0] * self.sigma_range[0])
        width_min = 1.0/(2.0 * self.sigma_range[1] * self.sigma_range[1])
        self.clib.daoInitialize2D(self.mfit,
                                  width_min,
                                  width_max)
        

class MultiFitter2DALS(MultiFitterGaussian):
    """
    Fit with a variable peak width (of the same size in X and Y) using
    the Anscombe least squares fitting model.
    """
    def __init__(self, sigma_range = None, **kwds):
        super(MultiFitter2DALS, self).__init__(**kwds)
        
        self.sigma_range = sigma_range
            
    def initializeC(self, image):
        super(MultiFitter2DALS, self).initializeC(image)

        width_max = 1.0/(2.0 * self.sigma_range[0] * self.sigma_range[0])
        width_min = 1.0/(2.0 * self.sigma_range[1] * self.sigma_range[1])
        self.clib.daoInitialize2DALS(self.mfit,
                                     width_min,
                                     width_max)

    def newImage(self, image):
         super(MultiFitter2DALS, self).newImage(image)
         self.clib.mFitAnscombeTransformImage(self.mfit)

         
class MultiFitter2DLS(MultiFitterGaussian):
    """
    Fit with a variable peak width (of the same size in X and Y) using the
    least squares fitting model.
    """
    def __init__(self, sigma_range = None, **kwds):
        super(MultiFitter2DLS, self).__init__(**kwds)
        
        self.sigma_range = sigma_range
            
    def initializeC(self, image):
        super(MultiFitter2DLS, self).initializeC(image)

        width_max = 1.0/(2.0 * self.sigma_range[0] * self.sigma_range[0])
        width_min = 1.0/(2.0 * self.sigma_range[1] * self.sigma_range[1])
        self.clib.daoInitialize2DLS(self.mfit,
                                    width_min,
                                    width_max)


class MultiFitter2DDWLS(MultiFitterGaussian):
    """
    Fit with a variable peak width (of the same size in X and Y) using the
    data weighted least squares fitting model.
    """
    def __init__(self, sigma_range = None, **kwds):
        super(MultiFitter2DDWLS, self).__init__(**kwds)
        
        self.sigma_range = sigma_range
            
    def initializeC(self, image):
        super(MultiFitter2DDWLS, self).initializeC(image)

        width_max = 1.0/(2.0 * self.sigma_range[0] * self.sigma_range[0])
        width_min = 1.0/(2.0 * self.sigma_range[1] * self.sigma_range[1])
        self.clib.daoInitialize2DDWLS(self.mfit,
                                      width_min,
                                      width_max)


class MultiFitter2DFWLS(MultiFitterGaussian):
    """
    Fit with a variable peak width (of the same size in X and Y) using the
    fit weighted least squares fitting model.
    """
    def __init__(self, sigma_range = None, **kwds):
        super(MultiFitter2DFWLS, self).__init__(**kwds)
        
        self.sigma_range = sigma_range
            
    def initializeC(self, image):
        super(MultiFitter2DFWLS, self).initializeC(image)

        width_max = 1.0/(2.0 * self.sigma_range[0] * self.sigma_range[0])
        width_min = 1.0/(2.0 * self.sigma_range[1] * self.sigma_range[1])
        self.clib.daoInitialize2DFWLS(self.mfit,
                                      width_min,
                                      width_max)
        

class MultiFitter3D(MultiFitterGaussian):
    """
    Fit with peak width that can change independently in X and Y.
    """
    def __init__(self, sigma_range = None, **kwds):
        super(MultiFitter3D, self).__init__(**kwds)
        
        self.sigma_range = sigma_range
        
    def initializeC(self, image):
        super(MultiFitter3D, self).initializeC(image)

        width_max = 1.0/(2.0 * self.sigma_range[0] * self.sigma_range[0])
        width_min = 1.0/(2.0 * self.sigma_range[1] * self.sigma_range[1])
        
        self.clib.daoInitialize3D(self.mfit,
                                  width_min,
                                  width_max)
        
        
class MultiFitterZ(MultiFitterGaussian):
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
# The MIT License
#
# Copyright (c) 2018 Zhuang Lab, Harvard University
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
