#!/usr/bin/env python
"""
Simple Python interface to cubic_fit.c.

Hazen 01/14
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.dao_fit_c as daoFitC

import storm_analysis.spliner.spline2D as spline2D
import storm_analysis.spliner.spline3D as spline3D


def loadCubicFitC():
    cubic_fit = loadclib.loadCLibrary("storm_analysis.spliner", "cubic_fit")

    # From sa_library/multi_fit.c
    cubic_fit.mFitGetFitImage.argtypes = [ctypes.c_void_p,
                                          ndpointer(dtype=numpy.float64)]

    cubic_fit.mFitGetNError.argtypes = [ctypes.c_void_p]
    cubic_fit.mFitGetNError.restype = ctypes.c_int
    
    cubic_fit.mFitGetPeakPropertyDouble.argtypes = [ctypes.c_void_p,
                                                    ndpointer(dtype=numpy.float64),
                                                    ctypes.c_char_p]
    
    cubic_fit.mFitGetPeakPropertyInt.argtypes = [ctypes.c_void_p,
                                                 ndpointer(dtype=numpy.int32),
                                                 ctypes.c_char_p]
    
    cubic_fit.mFitGetResidual.argtypes = [ctypes.c_void_p,
                                          ndpointer(dtype=numpy.float64)]

    cubic_fit.mFitGetUnconverged.argtypes = [ctypes.c_void_p]
    cubic_fit.mFitGetUnconverged.restype = ctypes.c_int

    cubic_fit.mFitIterateLM.argtypes = [ctypes.c_void_p]
    cubic_fit.mFitIterateOriginal.argtypes = [ctypes.c_void_p]

    cubic_fit.mFitNewBackground.argtypes = [ctypes.c_void_p,
                                            ndpointer(dtype=numpy.float64)]
    
    cubic_fit.mFitNewImage.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]

    cubic_fit.mFitRemoveErrorPeaks.argtypes = [ctypes.c_void_p]

    cubic_fit.mFitRemoveRunningPeaks.argtypes = [ctypes.c_void_p]

    cubic_fit.mFitSetPeakStatus.argtypes = [ctypes.c_void_p,
                                            ndpointer(dtype=numpy.int32)]    
    
    # From spliner/cubic_spline.c
    cubic_fit.getZSize.argtypes = [ctypes.c_void_p]
    cubic_fit.getZSize.restype = ctypes.c_int

    cubic_fit.initSpline2D.argtypes = [ndpointer(dtype=numpy.float64),
                                       ctypes.c_int,
                                       ctypes.c_int]
    cubic_fit.initSpline2D.restype = ctypes.c_void_p
    
    cubic_fit.initSpline3D.argtypes = [ndpointer(dtype=numpy.float64),
                                       ctypes.c_int,
                                       ctypes.c_int,
                                       ctypes.c_int]
    cubic_fit.initSpline3D.restype = ctypes.c_void_p

    # From spliner/cubic_fit.c
    cubic_fit.cfCleanup.argtypes = [ctypes.c_void_p]

    cubic_fit.cfInitialize.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64),
                                       ndpointer(dtype=numpy.float64),
                                       ctypes.c_double,
                                       ctypes.c_int,
                                       ctypes.c_int]
    cubic_fit.cfInitialize.restype = ctypes.POINTER(daoFitC.fitData)
    cubic_fit.cfInitialize2D.argtypes = [ctypes.c_void_p]
    cubic_fit.cfInitialize3D.argtypes = [ctypes.c_void_p]

    cubic_fit.cfNewPeaks.argtypes = [ctypes.c_void_p,
                                     ndpointer(dtype=numpy.float64),
                                     ctypes.c_char_p,
                                     ctypes.c_int]

    return cubic_fit
    

#
# Classes.
#

class CSplineFit(daoFitC.MultiFitterArbitraryPSF):

    def __init__(self, spline_fn = None, **kwds):
        super(CSplineFit, self).__init__(**kwds)
        self.spline_fn = spline_fn
        
        self.clib = loadCubicFitC()

    def cleanup(self, spacing = "  ", verbose = True):
        super(CSplineFit, self).cleanup(spacing = spacing,
                                        verbose = verbose)
        if self.mfit is not None:
            self.clib.cfCleanup(self.mfit)
            self.mfit = None
        self.c_spline = None

    def getSplineSize(self):
        return self.spline_fn.getSplineSize()
        
    def getSize(self):
        return self.spline_fn.getSize()
        
    def initializeC(self, image):
        """
        This initializes the C fitting library.
        """
        super(CSplineFit, self).initializeC(image)
        
        self.mfit = self.clib.cfInitialize(self.spline_fn.getCPointer(),
                                           self.scmos_cal,
                                           numpy.ascontiguousarray(self.clamp),
                                           self.default_tol,
                                           self.scmos_cal.shape[1],
                                           self.scmos_cal.shape[0])

    def newPeaks(self, peaks, peaks_type):
        """
        Pass new peaks to the C library.
        """
        c_peaks = self.formatPeaks(peaks, peaks_type)
        self.clib.cfNewPeaks(self.mfit,
                              c_peaks,
                              ctypes.c_char_p(peaks_type.encode()),
                              c_peaks.shape[0])


class CSpline2DFit(CSplineFit):

    def __init__(self, **kwds):
        super(CSpline2DFit, self).__init__(**kwds)
        
        # Clamp parameters.
        #
        # These are basically the same as the base class except for z.
        #
        self.clamp = numpy.array([1.0,  # Height (Note: This is relative to the initial guess).
                                  1.0,  # x position
                                  0.3,  # width in x
                                  1.0,  # y position
                                  0.3,  # width in y
                                  1.0,  # background (Note: This is relative to the initial guess).
                                  1.0]) # z position

    def initializeC(self, image):
        super(CSpline2DFit, self).initializeC(image)
        self.clib.cfInitialize2D(self.mfit)
        
    def rescaleZ(self, peaks):
        return peaks
        

class CSpline3DFit(CSplineFit):
    
    def __init__(self, **kwds):
        super(CSpline3DFit, self).__init__(**kwds)

        # Clamp parameters.
        #
        # Special clamp parameters for 3D. Basically we allow Z to change more, however
        # it is not obvious that this helps at all.
        #
        if True:
            self.clamp = numpy.array([1.0,  # Height (Note: This is relative to the initial guess).
                                      1.0,  # x position
                                      0.3,  # width in x
                                      1.0,  # y position
                                      0.3,  # width in y
                                      1.0,  # background (Note: This is relative to the initial guess).
                                      0.5 * self.spline_fn.getSplineSize()]) # z position (in spline size units).

    def initializeC(self, image):
        super(CSpline3DFit, self).initializeC(image)
        self.clib.cfInitialize3D(self.mfit)

    def rescaleZ(self, z):
        return self.spline_fn.rescaleZ(z)
