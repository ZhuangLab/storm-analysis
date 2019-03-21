#!/usr/bin/env python
"""
Simple Python interface to fft_fit.c.

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


def loadFFTFitC():
    fft_fit = loadclib.loadCLibrary("fft_fit")

    # From sa_library/multi_fit.c
    fft_fit.mFitGetFitImage.argtypes = [ctypes.c_void_p,
                                        ndpointer(dtype=numpy.float64)]


    fft_fit.mFitGetNError.argtypes = [ctypes.c_void_p]
    fft_fit.mFitGetNError.restype = ctypes.c_int
    
    fft_fit.mFitGetPeakPropertyDouble.argtypes = [ctypes.c_void_p,
                                                  ndpointer(dtype=numpy.float64),
                                                  ctypes.c_char_p]
    
    fft_fit.mFitGetPeakPropertyInt.argtypes = [ctypes.c_void_p,
                                               ndpointer(dtype=numpy.int32),
                                               ctypes.c_char_p]    

    fft_fit.mFitGetResidual.argtypes = [ctypes.c_void_p,
                                        ndpointer(dtype=numpy.float64)]

    fft_fit.mFitGetUnconverged.argtypes = [ctypes.c_void_p]
    fft_fit.mFitGetUnconverged.restype = ctypes.c_int

    fft_fit.mFitIterateLM.argtypes = [ctypes.c_void_p]

    fft_fit.mFitNewBackground.argtypes = [ctypes.c_void_p,
                                          ndpointer(dtype=numpy.float64)]    
    
    fft_fit.mFitNewImage.argtypes = [ctypes.c_void_p,
                                     ndpointer(dtype=numpy.float64)]

    fft_fit.mFitRemoveErrorPeaks.argtypes = [ctypes.c_void_p]

    fft_fit.mFitRemoveRunningPeaks.argtypes = [ctypes.c_void_p]

    fft_fit.mFitSetPeakStatus.argtypes = [ctypes.c_void_p,
                                          ndpointer(dtype=numpy.int32)]
    
    # From psf_fft/fft_fit.c
    fft_fit.ftFitCleanup.argtypes = [ctypes.c_void_p]

    fft_fit.ftFitInitialize.argtypes = [ctypes.c_void_p,
                                        ndpointer(dtype=numpy.float64),
                                        ndpointer(dtype=numpy.float64),
                                        ctypes.c_double,
                                        ctypes.c_int,
                                        ctypes.c_int]
    fft_fit.ftFitInitialize.restype = ctypes.POINTER(daoFitC.fitData)

    fft_fit.ftFitNewPeaks.argtypes = [ctypes.c_void_p,
                                      ndpointer(dtype=numpy.float64),
                                      ctypes.c_char_p,
                                      ctypes.c_int]

    return fft_fit
    

#
# Classes.
#
class CFFTFit(daoFitC.MultiFitterArbitraryPSF):

    def __init__(self, psf_fn = None, **kwds):
        super(CFFTFit, self).__init__(**kwds)
        self.psf_fn = psf_fn
        
        self.clib = loadFFTFitC()

    def cleanup(self, spacing = "  ", verbose = True):
        super(CFFTFit, self).cleanup(spacing = spacing,
                                     verbose = verbose)
        if self.mfit is not None:
            self.clib.ftFitCleanup(self.mfit)
            self.mfit = None
        self.psf_fft = None

    def getSize(self):
        return self.psf_fn.getSize()
        
    def initializeC(self, image):
        """
        This initializes the C fitting library.
        """
        super(CFFTFit, self).initializeC(image)

        self.mfit = self.clib.ftFitInitialize(self.psf_fn.getCPointer(),
                                              self.rqe,
                                              self.scmos_cal,
                                              self.default_tol,
                                              self.scmos_cal.shape[1],
                                              self.scmos_cal.shape[0])

    def newPeaks(self, peaks, peaks_type):
        """
        Pass new peaks to the C library.
        """
        c_peaks = self.formatPeaks(peaks, peaks_type)
        self.clib.ftFitNewPeaks(self.mfit,
                                c_peaks,
                                ctypes.c_char_p(peaks_type.encode()),
                                c_peaks.shape[0])

    def rescaleZ(self, z):
        return self.psf_fn.rescaleZ(z)


