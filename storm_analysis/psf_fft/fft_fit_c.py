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
    fft_fit = loadclib.loadCLibrary("storm_analysis.pupilfn", "fft_fit")

    # From sa_library/multi_fit.c
    fft_fit.mFitGetFitImage.argtypes = [ctypes.c_void_p,
                                        ndpointer(dtype=numpy.float64)]

    fft_fit.mFitGetResidual.argtypes = [ctypes.c_void_p,
                                        ndpointer(dtype=numpy.float64)]
    
    fft_fit.mFitGetResults.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]

    fft_fit.mFitGetUnconverged.argtypes = [ctypes.c_void_p]
    fft_fit.mFitGetUnconverged.restype = ctypes.c_int

    fft_fit.mFitIterateLM.argtypes = [ctypes.c_void_p]
    fft_fit.mFitIterateOriginal.argtypes = [ctypes.c_void_p]
    
    fft_fit.mFitNewImage.argtypes = [ctypes.c_void_p,
                                     ndpointer(dtype=numpy.float64)]
    
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
                                      ctypes.c_int]

    return fft_fit
    

#
# Classes.
#
class CFFTFit(daoFitC.MultiFitterBase):

    def __init__(self, psf_fn = None, **kwds):
        super(FFTFit, self).__init__(**kwds)
        
        self.psf_fn = psf_fn

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
        super(CFFTFit, self).cleanup(spacing = spacing,
                                     verbose = verbose)
        if self.mfit is not None:
            self.clib.ftFitCleanup(self.mfit)
            self.mfit = None
        self.psf_fft = None

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
        return self.psf_fn.getSize()
        
    def initializeC(self, image):
        """
        This initializes the C fitting library. You can call this directly, but
        the idea is that it will get called automatically the first time that you
        provide a new image for fitting.
        """
        super(CFFTFit, self).initializeC(image)

        self.mfit = self.clib.ftFitInitialize(self.psf_fn.getCPointer(),
                                              self.scmos_cal,
                                              numpy.ascontiguousarray(self.clamp),
                                              self.default_tol,
                                              self.scmos_cal.shape[1],
                                              self.scmos_cal.shape[0])

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

    def rescaleZ(self, peaks):
        z_index = utilC.getZCenterIndex()
        peaks[:,z_index] = self.psf_fn.rescaleZ(peaks[:,z_index])
        return peaks


