#!/usr/bin/env python
"""
Python interface to mp_fit_arb.c. 

Hazen 01/18
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer

import storm_analysis.sa_library.dao_fit_c as daoFitC
import storm_analysis.sa_library.loadclib as loadclib

import storm_analysis.multi_plane.mp_fit_c as mpFitC


def loadMPFitCArb():
    mp_fit = loadclib.loadCLibrary("mp_fit_arb")
    
    # From multi_plane/mp_fit_arb.c
    mp_fit.mpArbInitializePSFFFTChannel.argtypes = [ctypes.c_void_p,
                                                    ctypes.c_void_p,
                                                    ndpointer(dtype=numpy.float64),
                                                    ndpointer(dtype=numpy.float64),
                                                    ctypes.c_int]
    
    mp_fit.mpArbInitializePupilFnChannel.argtypes = [ctypes.c_void_p,
                                                     ctypes.c_void_p,
                                                     ndpointer(dtype=numpy.float64),
                                                     ndpointer(dtype=numpy.float64),
                                                     ctypes.c_double,
                                                     ctypes.c_double,
                                                     ctypes.c_int]
    
    mp_fit.mpArbInitializeSplineChannel.argtypes = [ctypes.c_void_p,
                                                    ctypes.c_void_p,
                                                    ndpointer(dtype=numpy.float64),
                                                    ndpointer(dtype=numpy.float64),
                                                    ctypes.c_int]

    mp_fit.mpArbNewPeaks.argtypes = [ctypes.c_void_p,
                                     ndpointer(dtype=numpy.float64),
                                     ctypes.c_char_p,
                                     ctypes.c_int]    
    
    return mp_fit


class MPFitArb(mpFitC.MPFit):
    """
    Base class for multi-plane fitting with an arbitrary PSF.
    """
    def __init__(self, independent_heights = None, psf_objects = None, **kwds):
        super(MPFitArb, self).__init__(**kwds)

        self.independent_heights = independent_heights
        self.n_channels = len(psf_objects)
        self.psf_objects = psf_objects

        self.clib = loadMPFitCArb()
        mpFitC.addMPFitC(self.clib)

    def getSize(self):
        return self.psf_objects[0].getSize()        

    def newPeaks(self, peaks, peaks_type):
        c_peaks = daoFitC.formatPeaksArbitraryPSF(peaks, peaks_type)
        self.clib.mpArbNewPeaks(self.mfit,
                                c_peaks,
                                ctypes.c_char_p(peaks_type.encode()),
                                c_peaks.shape[0])


class MPPSFFnFit(MPFitArb):
    """
    The basic idea is that we are going to use the functionality from PSF FFT
    to do most of the work. We will have one psfFFTFit C structure per image
    plane / channel.
    """
    def __init__(self, **kwds):
        super(MPPSFFnFit, self).__init__(**kwds)
    
    def initializeChannel(self, rqe, variance, channel):
        super(MPPSFFnFit, self).initializeChannel(rqe, variance, channel)
        
        # This where the differentation in which type of fitter to use happens.
        zmax = self.psf_objects[0].getZMax() * 1.0e-3
        zmin = self.psf_objects[0].getZMin() * 1.0e-3
        self.clib.mpArbInitializePSFFFTChannel(self.mfit,
                                               self.psf_objects[channel].getCPointer(),
                                               rqe,
                                               variance,
                                               channel)

    def rescaleZ(self, z):
        return self.psf_objects[0].rescaleZ(z)
        
    def setWeights(self, weights, verbose = True):
        if weights is None:
            weights = {"bg" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "h" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "x" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "y" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "z" : numpy.ones((2, self.n_channels))/float(self.n_channels)}
            super(MPSFFnFit, self).setWeights(weights, 0.0, 0.0, verbose = verbose)

        else:
            zmax = self.psf_objects[0].getZMax() * 1.0e-3
            zmin = self.psf_objects[0].getZMin() * 1.0e-3
            z_offset = -0.5*float(self.psf_objects[0].getZSize())
            z_scale = float(weights["bg"].shape[0]-1)/float(self.psf_objects[0].getZSize())
            super(MPPSFFnFit, self).setWeights(weights, z_offset, z_scale, verbose = verbose)

            
class MPPupilFnFit(MPFitArb):
    """
    The basic idea is that we are going to use the functionality from PupilFn
    to do most of the work. We will have one pupilFit C structure per image
    plane / channel.
    """
    def __init__(self, **kwds):
        super(MPPupilFnFit, self).__init__(**kwds)

    def initializeChannel(self, rqe, variance, channel):
        super(MPPupilFnFit, self).initializeChannel(rqe, variance, channel)
        
        # This where the differentation in which type of fitter to use happens.
        zmax = self.psf_objects[0].getZMax() * 1.0e-3
        zmin = self.psf_objects[0].getZMin() * 1.0e-3
        self.clib.mpArbInitializePupilFnChannel(self.mfit,
                                                self.psf_objects[channel].getCPointer(),
                                                rqe,
                                                variance,
                                                zmin,
                                                zmax,
                                                channel)

    def setWeights(self, weights):
        if weights is None:
            weights = {"bg" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "h" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "x" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "y" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "z" : numpy.ones((2, self.n_channels))/float(self.n_channels)}
            super(MPPupilFnFit, self).setWeights(weights, 0.0, 0.0)

        else:
            zmax = self.psf_objects[0].getZMax() * 1.0e-3
            zmin = self.psf_objects[0].getZMin() * 1.0e-3
            z_offset = zmin
            z_scale = float(weights["bg"].shape[0]-1)/(zmax - zmin + 1.0e-12)
            super(MPPupilFnFit, self).setWeights(weights, z_offset, z_scale)
            
    
class MPSplineFit(MPFitArb):
    """
    The basic idea is that we are going to use the functionality from Spliner
    to do most of the work. We will have one splineFit C structure per image
    plane / channel.
    """
    def __init__(self, **kwds):
        super(MPSplineFit, self).__init__(**kwds)
    
    def initializeChannel(self, rqe, variance, channel):
        super(MPSplineFit, self).initializeChannel(rqe, variance, channel)
        
        # This where the differentation in which type of fitter to use happens.
        self.clib.mpArbInitializeSplineChannel(self.mfit,
                                               self.psf_objects[channel].getCPointer(),
                                               rqe,
                                               variance,
                                               channel)

    def rescaleZ(self, z):
        return self.psf_objects[0].rescaleZ(z)
        
    def setWeights(self, weights):
        if weights is None:
            weights = {"bg" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "h" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "x" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "y" : numpy.ones((2, self.n_channels))/float(self.n_channels),
                       "z" : numpy.ones((2, self.n_channels))/float(self.n_channels)}
            super(MPSplineFit, self).setWeights(weights, 0.0, 0.0)

        else:
            assert (weights["bg"].shape[0] == (self.psf_objects[0].getSize() + 1)), "Incorrect shape for weights array."
            super(MPSplineFit, self).setWeights(weights, 0.0, 1.0)
            
