#!/usr/bin/env python
"""
Simple Python interface to mp_fit.c.

Hazen 01/14
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.dao_fit_c as daoFitC
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.loadclib as loadclib

import storm_analysis.spliner.spline3D as spline3D

def loadMPFitC():
    mp_fit = loadclib.loadCLibrary("storm_analysis.multi_plane", "mp_fit")

    # From spliner/cubic_spline.c
    mp_fit.getZSize.argtypes = [ctypes.c_void_p]
    mp_fit.getZSize.restype = ctypes.c_int
    
    mp_fit.initSpline3D.argtypes = [ndpointer(dtype=numpy.float64),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]
    mp_fit.initSpline3D.restype = ctypes.c_void_p

    # From multi_plane/mp_fit.c
    mp_fit.mpCleanup.argtypes = [ctypes.c_void_p]

    mp_fit.mpGetFitImage.argtypes = [ctypes.c_void_p,
                                     ndpointer(dtype=numpy.float64),
                                     ctypes.c_int]

    mp_fit.mpGetResults.argtypes = [ctypes.c_void_p,
                                    ndpointer(dtype=numpy.float64)]

    mp_fit.mpGetUnconverged.argtypes = [ctypes.c_void_p]
    mp_fit.mpGetUnconverged.restype = ctypes.c_int
    
    mp_fit.mpInitialize.argtypes = [ndpointer(dtype=numpy.float64),
                                    ctypes.c_double,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]
    mp_fit.mpInitialize.restype = ctypes.c_void_p

    mp_fit.mpInitializeChannel.argtypes = [ctypes.c_void_p,
                                           ctypes.c_void_p,
                                           ndpointer(dtype=numpy.float64),
                                           ctypes.c_int]
    
    mp_fit.mpIterate.argtypes = [ctypes.c_void_p]

    mp_fit.mpNewImage.argtypes = [ctypes.c_void_p,
                                  ndpointer(dtype=numpy.float64),
                                  ctypes.c_int]

    mp_fit.mpNewPeaks.argtypes = [ctypes.c_void_p,
                                  ndpointer(dtype=numpy.float64),
                                  ctypes.c_int]

    mp_fit.mpSetTransforms.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64),
                                       ndpointer(dtype=numpy.float64),
                                       ndpointer(dtype=numpy.float64),
                                       ndpointer(dtype=numpy.float64)]

    return mp_fit


class MPSplineFit(daoFitC.MultiFitterBase):
    """
    Multi-plane fitting object. While technically a sub-class of MultiFitterBase
    this overrides essentially all of the base class methods.

    The basic idea is that we are going to use the functionality from Spliner
    to do most of the work. We will have one splineFit C structure per image
    plane / channel.
    """
    def __init__(self, splines, coeffs, verbose = False):
        super().__init__(verbose)

        self.clib = loadMPFitC()
        self.c_splines = []
        self.mappings = None
        self.n_channels = 0
        self.py_splines = []
        self.xt = []
        self.yt = []

        # Default clamp parameters.
        #
        # These set the (initial) scale for how much these parameters
        # can change in a single fitting iteration.
        #
        self.clamp = numpy.array([1000.0,   # Height
                                  1.0,      # x position
                                  0.3,      # width in x
                                  1.0,      # y position
                                  0.3,      # width in y
                                  100.0,    # background
                                  1.0])     # z position

        # Initialize splines.
        for i in range(len(splines)):
            self.py_splines.append(spline3D.Spline3D(splines[i], coeff = coeffs[i]))
            self.c_splines.append(self.clib.initSpline3D(numpy.ascontiguousarray(self.py_splines[i].coeff, dtype = numpy.float64),
                                                         self.py_splines[i].max_i,
                                                         self.py_splines[i].max_i,
                                                         self.py_splines[i].max_i))
            self.n_channels += 1
        self.inv_zscale = 1.0/self.clib.getZSize(self.c_splines[0])

        # Initialize storage for mappings.
        for i in range(2):
            self.xt.append(numpy.zeros((self.n_channels, 3)))
            self.yt.append(numpy.zeros((self.n_channels, 3)))

    def cleanup(self):
        if self.mfit is not None:
            self.clib.cleanup(self.mfit)

    def doFit(self, peaks, max_iterations = 200):

        # Initialize C library with new peaks.
        self.clib.mpNewPeaks(self.mfit,
                             numpy.ascontiguousarray(peaks),
                             peaks.shape[0])

        # Iterate fittings.
        i = 0
        self.iterate()
        while(self.clib.mpGetUnconverged(self.mfit) and (i < max_iterations)):
            if self.verbose and ((i%20)==0):
                print("iteration", i)
            self.iterate()
            i += 1

        if self.verbose:
            if (i == max_iterations):
                print(" Failed to converge in:", i, self.clib.mpGetUnconverged(self.mfit))
            else:
                print(" Multi-fit converged in:", i, self.clib.mpGetUnconverged(self.mfit))
            print("")

        # Get updated peak values back from the C library.
        fit_peaks = numpy.ascontiguousarray(numpy.zeros(peaks.shape))
        self.clib.mpGetResults(self.mfit, fit_peaks)

        return fit_peaks            

    def getFitImage(self, channel):
        fit_image = numpy.zeros(self.im_shape)
        self.clib.mpGetFitImage(self.mfit, fit_image, channel)
        return fit_image

    def getGoodPeaks(self, peaks, threshold):
        if (peaks.size > 0):
            #
            # The peaks are in groups that is n_channels long, one peak
            # per channel. In the C library, all peaks in the same group
            # will have the same status index and height, so we can filter
            # here with some confidence that we'll keep or eliminate all
            # the peaks in a particular group at the same time. We attempt
            # to enforce this with an assertion statement, but it is an
            # imperfect test as we could eliminate 1 peak from each of 4
            # groups and still pass, even though this would completely
            # mess up the analysis.
            #
            status_index = utilC.getStatusIndex()
            height_index = utilC.getHeightIndex()

            mask = (peaks[:,status_index] != 2.0) & (peaks[:,height_index] > min_height)
            if self.verbose:
                print(" ", numpy.sum(mask), "were good out of", peaks.shape[0])

            masked_peaks = peaks[mask,:]

            # The number of peaks should always be a multiple of the number of channels.
            assert((masked_peaks.shape[0] % self.n_channels) == 0)
                
            return masked_peaks
                
        else:
            return peaks
    
    def initializeC(self, variance):
        self.im_shape = variance.shape

        # Get fitting structure.
        self.mfit = self.clib.mpInitialize(numpy.ascontiguousarray(self.clamp),
                                           self.default_tol,
                                           self.n_channels,
                                           variance.shape[1],
                                           variance.shape[0])

        # Add transforms.
        self.clib.mpSetTransforms(self.mfit,
                                  numpy.ascontiguousarray(self.xt[0]),
                                  numpy.ascontiguousarray(self.yt[0]),
                                  numpy.ascontiguousarray(self.xt[1]),
                                  numpy.ascontiguousarray(self.yt[1]))
    
    def iterate(self):
        self.clib.mpIterate(self.mfit)

    def newImage(self, image, channel):
        if (image.shape[0] != self.im_shape[0]) or (image.shape[1] != self.im_shape[1]):
            raise daoFitC.MultiFitterException("Current image shape and the original image shape are not the same.")

        self.clib.mpNewImage(self.mfit, image, channel)

    def rescaleZ(self, peaks, zmin, zmax):
        z_index = utilC.getZCenterIndex()
        spline_range = zmax - zmin
        peaks[:,z_index] = peaks[:,z_index] * self.inv_zscale * spline_range + zmin
        return peaks

    def setMapping(self, ch_from, ch_to, xt, yt):
        #
        # We temporarily store the channel to channel affine transform
        # mapping coefficients until we have initialized the C library.
        #

        # These are the transforms to go from channel 0 to channel N.
        if (ch_from == 0):
            self.xt[0][ch_to,:] = xt
            self.yt[0][ch_to,:] = yt

        # These are the transforms to go from channel N to channel 0.
        else:
            self.xt[1][ch_from,:] = xt
            self.yt[1][ch_from,:] = yt

    def setVariance(self, variance, channel):
        #
        # This is a little difference than the 3D-DAOSTORM, sCMOS and Spliner
        # because this is the first thing that will be called with an array
        # that is expected to have the same size of the images. So we initialize
        # the C library now, rather than in newImage().
        #
        if self.mfit is None:
            self.initializeC(variance)
        else:
            if (variance.shape[0] != self.im_shape[0]) or (variance.shape[1] != self.im_shape[1]):
                raise daoFitC.MultiFitterException("Current variance shape and the original variance shape are not the same.")

        self.clib.mpInitializeChannel(self.mfit,
                                      self.c_splines[channel],
                                      variance,
                                      channel)


        
