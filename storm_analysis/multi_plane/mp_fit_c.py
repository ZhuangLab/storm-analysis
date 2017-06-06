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

    # C interface definition.

    # From sa_library/multi_fit.c
    mp_fit.mFitGetResults.argtypes = [ctypes.c_void_p,
                                      ndpointer(dtype=numpy.float64)]

    mp_fit.mFitGetUnconverged.argtypes = [ctypes.c_void_p]
    mp_fit.mFitGetUnconverged.restype = ctypes.c_int

    # From spliner/cubic_spline.c
    mp_fit.initSpline3D.argtypes = [ndpointer(dtype=numpy.float64),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]
    mp_fit.initSpline3D.restype = ctypes.c_void_p

    # From multi_plane/mp_fit.c
    mp_fit.cleanup.argtypes = [ctypes.c_void_p]

    mp_fit.getFitImage.argtypes = [ctypes.c_void_p,
                                   ndpointer(dtype=numpy.float64),
                                   ctypes.c_int]
    
    mp_fit.initialize.argtypes = [ndpointer(dtype=numpy.float64),
                                  ctypes.c_double,
                                  ctypes.c_int,
                                  ctypes.c_int]
    mp_fit.initialize.restype = ctypes.c_void_p

    mp_fit.iterate.argtypes = [ctypes.c_void_p]

    mp_fit.newImage.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype=numpy.float64),
                                ctypes.c_int]

    mp_fit.newPeaks.argtypes = [ctypes.c_void_p,
                                ndpointer(dtype=numpy.float64),
                                ctypes.c_int]

    mp_fit.setCSpline.argtypes = [ctypes.c_void_p,
                                  ctypes.c_void_p,
                                  ctypes.c_int]

    mp_fit.setTransforms.argtypes = [ctypes.c_void_p,
                                     ndpointer(dtype=numpy.float64),
                                     ndpointer(dtype=numpy.float64),
                                     ndpointer(dtype=numpy.float64),
                                     ndpointer(dtype=numpy.float64)]

    mp_fit.setVariance.argtypes = [ctypes.c_void_p,
                                   ndpointer(dtype=numpy.float64),
                                   ctypes.c_int]    

    return mp_fit


class MPSplineFit(daoFitC.MultiFitterBase):

    def __init__(self, splines, coeffs, verbose = False):
        super().__init__(verbose)

        self.clib = loadMPFitC()
        self.c_splines = []
        self.mappings = None
        self.n_channels = 0
        self.py_splines = []
        self.xt = []
        self.yt = []

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
            self.xt.append(numpy.zeros((self.n_channels-1, 3)))
            self.yt.append(numpy.zeros((self.n_channels-1, 3)))

    def cleanup(self):
        super().cleanup()

        #
        # The C library will free the spline storage, we just need to
        # reset the pointers so that we can't use them by accident. Not
        # that we would be able to anyway..
        #
        for i in range(len(self.c_splines))
            self.c_splines = None
        
    def getFitImage(self, channel):
        fit_image = numpy.zeros(self.im_shape)
        self.clib.getFitImage(self.mfit, fit_image, channel)
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
        self.mfit = self.clib.initialize(numpy.ascontiguousarray(self.clamp),
                                         self.default_tol,
                                         variance.shape[1],
                                         variance.shape[0])

        # Add splines.
        for i in range(len(self.c_splines)):
            self.clib.setCSpline(self.mfit, self.c_splines[i], i)

        # Add transforms.
        self.clib.setTransforms(self.mfit,
                                numpy.ascontiguousarray(self.xt[0]),
                                numpy.ascontiguousarray(self.yt[0]),
                                numpy.ascontiguousarray(self.xt[1]),
                                numpy.ascontiguousarray(self.yt[1]))
    
    def iterate(self):
        self.clib.iterate(self.mfit)

    def newImage(self, image, channel):
        if (image.shape[0] != self.im_shape[0]) or (image.shape[1] != self.im_shape[1]):
            raise daoFitC.MultiFitterException("Current image shape and the original image shape are not the same.")

        self.clib.newImage(self.mfit, image, channel)

    def rescaleZ(self, peaks, zmin, zmax):
        z_index = utilC.getZCenterIndex()
        spline_range = zmax - zmin
        peaks[:,z_index] = peaks[:,z_index] * self.inv_zscale * spline_range + zmin
        return peaks

    def setMapping(self, ch_from, ch_to, xt, yt):
        #
        # We temporarily store the channel to channel affine transform
        # mapping coefficients until we have initialize the C library.
        #

        # These are the transforms to go from channel 0 to some other channel.
        if (ch_from == 0):
            self.xt[0][(ch_to-1),:] = xt
            self.yt[0][(ch_to-1),:] = yt

        # These are the transforms to go from some other channel to channel 0.
        else:
            self.xt[1][(ch_from-1),:] = xt
            self.yt[1][(ch_from-1),:] = yt

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

        self.clib.setVariance(self.mfit, variance, channel)


        
