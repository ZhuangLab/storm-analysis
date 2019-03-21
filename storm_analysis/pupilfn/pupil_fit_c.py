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
    pupil_fit = loadclib.loadCLibrary("pupil_fit")

    # From sa_library/multi_fit.c
    pupil_fit.mFitGetFitImage.argtypes = [ctypes.c_void_p,
                                          ndpointer(dtype=numpy.float64)]

    pupil_fit.mFitGetNError.argtypes = [ctypes.c_void_p]
    pupil_fit.mFitGetNError.restype = ctypes.c_int
    
    pupil_fit.mFitGetPeakPropertyDouble.argtypes = [ctypes.c_void_p,
                                                    ndpointer(dtype=numpy.float64),
                                                    ctypes.c_char_p]
    
    pupil_fit.mFitGetPeakPropertyInt.argtypes = [ctypes.c_void_p,
                                                 ndpointer(dtype=numpy.int32),
                                                 ctypes.c_char_p]

    pupil_fit.mFitGetResidual.argtypes = [ctypes.c_void_p,
                                          ndpointer(dtype=numpy.float64)]
    
    pupil_fit.mFitGetUnconverged.argtypes = [ctypes.c_void_p]
    pupil_fit.mFitGetUnconverged.restype = ctypes.c_int

    pupil_fit.mFitIterateLM.argtypes = [ctypes.c_void_p]

    pupil_fit.mFitNewBackground.argtypes = [ctypes.c_void_p,
                                            ndpointer(dtype=numpy.float64)]    
    
    pupil_fit.mFitNewImage.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]

    pupil_fit.mFitRemoveErrorPeaks.argtypes = [ctypes.c_void_p]

    pupil_fit.mFitRemoveRunningPeaks.argtypes = [ctypes.c_void_p]

    pupil_fit.mFitSetPeakStatus.argtypes = [ctypes.c_void_p,
                                            ndpointer(dtype=numpy.int32)]       
    
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
                                       ctypes.c_char_p,
                                       ctypes.c_int]

    pupil_fit.pfitSetZRange.argtypes = [ctypes.c_void_p,
                                        ctypes.c_double,
                                        ctypes.c_double]

    return pupil_fit
    

#
# Classes.
#
class CPupilFit(daoFitC.MultiFitterArbitraryPSF):

    def __init__(self, pupil_fn = None, **kwds):
        super(CPupilFit, self).__init__(**kwds)
        self.pupil_fn = pupil_fn
        
        self.clib = loadPupilFitC()

    def cleanup(self, spacing = "  ", verbose = True):
        super(CPupilFit, self).cleanup(spacing = spacing,
                                       verbose = verbose)
        if self.mfit is not None:
            self.clib.pfitCleanup(self.mfit)
            self.mfit = None
        self.pupil_fn = None

    def getSize(self):
        return self.pupil_fn.getSize()
        
    def initializeC(self, image):
        """
        This initializes the C fitting library.
        """
        super(CPupilFit, self).initializeC(image)

        self.mfit = self.clib.pfitInitialize(self.pupil_fn.getCPointer(),
                                             self.rqe,
                                             self.scmos_cal,
                                             self.default_tol,
                                             self.scmos_cal.shape[1],
                                             self.scmos_cal.shape[0])
        self.clib.pfitSetZRange(self.mfit,
                                self.pupil_fn.getZMin(),
                                self.pupil_fn.getZMax())

    def newPeaks(self, peaks, peaks_type):
        """
        Pass new peaks to the C library.
        """
        c_peaks = self.formatPeaks(peaks, peaks_type)
        self.clib.pfitNewPeaks(self.mfit,
                               c_peaks,
                               ctypes.c_char_p(peaks_type.encode()),
                               c_peaks.shape[0])
