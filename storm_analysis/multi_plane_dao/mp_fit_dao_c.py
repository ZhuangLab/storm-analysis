#!/usr/bin/env python
"""
Python interface to mp_fit_dao.c. 

Hazen 01/18
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer

import storm_analysis.sa_library.dao_fit_c as daoFitC
import storm_analysis.sa_library.loadclib as loadclib

import storm_analysis.multi_plane.mp_fit_c as mpFitC


def loadMPFitCDao():
    mp_fit = loadclib.loadCLibrary("storm_analysis.multi_plane", "mp_fit_dao")
    
    # From multi_plane/mp_fit_dao.c
    mp_fit.mpDaoInitialize2DChannel.argtypes = [ctypes.c_void_p,
                                                ndpointer(dtype=numpy.float64),
                                                ctypes.c_int,
                                                ctypes.c_int]

    mp_fit.mpDaoNewPeaks.argtypes = [ctypes.c_void_p,
                                     ndpointer(dtype=numpy.float64),
                                     ctypes.c_char_p,
                                     ctypes.c_int]    
    
    return mp_fit


class MPFitDao(mpFitC.MPFit):
    """
    Base class for multi-plane fitting.
    """
    def __init__(self, n_channels = None, roi_size = None, **kwds):
        super(MPFitDao, self).__init__(**kwds)

        self.independent_heights = 1
        self.n_channels = n_channels
        self.roi_size = roi_size

        self.clib = loadMPFitCDao()
        mpFitC.addMPFitC(self.clib)

    def newPeaks(self, peaks, peaks_type):
        c_peaks = daoFitC.formatPeaksGaussianPSF(peaks, peaks_type)
        self.clib.mpDaoNewPeaks(self.mfit,
                                c_peaks,
                                ctypes.c_char_p(peaks_type.encode()),
                                c_peaks.shape[0])


class MPFitDao2D(MPFitDao):
    """
    The basic idea is that we are going to use the functionality from 3D-DAOSTORM
    to do most of the work. We will have one 2D Gaussian fitter per image plane / 
    channel.
    """
    def __init__(self, **kwds):
        super(MPFitDao2D, self).__init__(**kwds)
    
    def setVariance(self, variance, channel):
        super(MPFitDao2D, self).setVariance(variance, channel)
        self.clib.mpDaoInitialize2DChannel(self.mfit,
                                           variance,
                                           self.roi_size,
                                           channel)

    def setWeights(self):
        """
        We can't weight by Z since we are not directly fitting for it, so use
        equal values for all the weights.
        """
        weights = {"bg" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                   "h" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                   "x" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                   "y" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                   "z" : numpy.ones((2, self.n_channels))/float(self.n_channels)}
        super(MPFitDao2D, self).setWeights(weights, 0.0, 0.0)
