#!/usr/bin/env python
"""
Simple Python interface to pupil_fit.c.

Hazen 10/17
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.dao_fit_c as daoFitC


def loadPupilFitC():
    pupil_fit = loadclib.loadCLibrary("storm_analysis.pupilfn", "pupil_fit")

    # From sa_library/multi_fit.c
    pupil_fit.mFitGetFitImage.argtypes = [ctypes.c_void_p,
                                          ndpointer(dtype=numpy.float64)]

    pupil_fit.mFitGetResidual.argtypes = [ctypes.c_void_p,
                                          ndpointer(dtype=numpy.float64)]
    
    pupil_fit.mFitGetResults.argtypes = [ctypes.c_void_p,
                                         ndpointer(dtype=numpy.float64)]

    pupil_fit.mFitGetUnconverged.argtypes = [ctypes.c_void_p]
    pupil_fit.mFitGetUnconverged.restype = ctypes.c_int

    pupil_fit.mFitIterateLM.argtypes = [ctypes.c_void_p]
    pupil_fit.mFitIterateOriginal.argtypes = [ctypes.c_void_p]
    
    pupil_fit.mFitNewImage.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]
    
    # From pupilfn/pupil_fit.c
    pupil_fit.pfitCleanup.argtypes = [ctypes.c_void_p]

    pupil_fit.pfitInitialize.argtypes = [ctypes.c_void_p,
                                         ndpointer(dtype=numpy.float64),
                                         ndpointer(dtype=numpy.float64),
                                         ctypes.c_double,
                                         ctypes.c_int,
                                         ctypes.c_int]
    pupil_fit.pfitInitialize.restype = ctypes.POINTER(daoFitC.fitData)

    pupil_fit.pfitNewPeaks.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64),
                                       ctypes.c_int]

    pupil_fit.pfitSetZRange.argtypes = [ctypes.c_void_p,
                                        ctypes.c_double,
                                        ctypes.c_double]

    return pupil_fit
    

#
# Classes.
#
class CPupilFit(daoFitC.MultiFitterBase):

    def __init__(self, pupil_fn = None, **kwds):
        super(CPupilFit, self).__init__(**kwds)
        self.pupil_fn = pupil_fn

        # Default clamp parameters.
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
        
        self.clib = loadPupilFitC()

    def cleanup(self, spacing = "  ", verbose = True):
        super(CPupilFit, self).cleanup(spacing = spacing,
                                       verbose = verbose)
        if self.mfit is not None:
            self.clib.pfitCleanup(self.mfit)
            self.mfit = None
        self.pupil_fn = None

    def getGoodPeaks(self, peaks, min_width):
        """
        min_width is ignored, it only exists so that this function has the correct signature.
        """
        if (peaks.size > 0):
            status_index = utilC.getStatusIndex()

            mask = (peaks[:,status_index] != 2.0)
            if self.verbose:
                print(" ", numpy.sum(mask), "were good out of", peaks.shape[0])
            return peaks[mask,:]
        else:
            return peaks

    def getSize(self):
        return self.pupil_fn.getSize()
        
    def initializeC(self, image):
        """
        This initializes the C fitting library. You can call this directly, but
        the idea is that it will get called automatically the first time that you
        provide a new image for fitting.
        """
        super(CPupilFit, self).initializeC(image)

        self.mfit = self.clib.pfitInitialize(self.pupil_fn.getCPointer(),
                                             self.scmos_cal,
                                             numpy.ascontiguousarray(self.clamp),
                                             self.default_tol,
                                             self.scmos_cal.shape[1],
                                             self.scmos_cal.shape[0])
        self.clib.pfitSetZRange(self.mfit, self.min_z, self.max_z)

    def iterate(self):
        self.clib.mFitIterateLM(self.mfit)
        #self.clib.mFitIterateOriginal(self.mfit)

    def newPeaks(self, peaks):
        """
        Pass new peaks to the C library.
        """
        self.clib.pfitNewPeaks(self.mfit,
                               numpy.ascontiguousarray(peaks),
                               peaks.shape[0])


