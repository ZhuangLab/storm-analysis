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

class mpFitData(ctypes.Structure):
    _fields_ = [('im_size_x', ctypes.c_int),
                ('im_size_y', ctypes.c_int),

                ('n_channels', ctypes.c_int),
                ('n_weights', ctypes.c_int),
                
                ('nfit', ctypes.c_int),
                
                ('w_z_offset', ctypes.c_double),
                ('w_z_scale', ctypes.c_double),
                
                ('tolerance', ctypes.c_double),
                
                ('clamp_start', (ctypes.c_double*7)),

                ('xt_0toN', ctypes.POINTER(ctypes.c_double)),
                ('yt_0toN', ctypes.POINTER(ctypes.c_double)),
                ('xt_Nto0', ctypes.POINTER(ctypes.c_double)),
                ('yt_Nto0', ctypes.POINTER(ctypes.c_double)),

                ('w_bg', ctypes.POINTER(ctypes.c_double)),
                ('w_h', ctypes.POINTER(ctypes.c_double)),
                ('w_x', ctypes.POINTER(ctypes.c_double)),
                ('w_y', ctypes.POINTER(ctypes.c_double)),
                ('w_z', ctypes.POINTER(ctypes.c_double)),
                ('heights', ctypes.POINTER(ctypes.c_double)),

                ('jacobian', ctypes.c_void_p),
                ('w_jacobian', ctypes.c_void_p),
                ('hessian', ctypes.c_void_p),
                ('w_hessian', ctypes.c_void_p),
                
                ('fit_data', ctypes.POINTER(ctypes.POINTER(daoFitC.fitData))),
                ('fn_update', ctypes.c_void_p)]

    
def loadMPFitC():
    mp_fit = loadclib.loadCLibrary("storm_analysis.multi_plane", "mp_fit")

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
                                    ctypes.c_int,
                                    ctypes.c_int]
    mp_fit.mpInitialize.restype = ctypes.POINTER(mpFitData)
#    mp_fit.mpInitialize.restype = ctypes.c_void_p

    mp_fit.mpInitializePupilFnChannel.argtypes = [ctypes.c_void_p,
                                                  ctypes.c_void_p,
                                                  ndpointer(dtype=numpy.float64),
                                                  ctypes.c_double,
                                                  ctypes.c_double,
                                                  ctypes.c_int]
    
    mp_fit.mpInitializeSplineChannel.argtypes = [ctypes.c_void_p,
                                                 ctypes.c_void_p,
                                                 ndpointer(dtype=numpy.float64),
                                                 ctypes.c_int]
    
    mp_fit.mpIterateLM.argtypes = [ctypes.c_void_p]
    mp_fit.mpIterateOriginal.argtypes = [ctypes.c_void_p]

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

    mp_fit.mpSetWeights.argtypes = [ctypes.c_void_p,
                                    ndpointer(dtype=numpy.float64),
                                    ndpointer(dtype=numpy.float64),
                                    ndpointer(dtype=numpy.float64),
                                    ndpointer(dtype=numpy.float64),
                                    ndpointer(dtype=numpy.float64),
                                    ctypes.c_int]
    
    mp_fit.mpSetWeightsIndexing.argtypes = [ctypes.c_void_p,
                                            ctypes.c_double,
                                            ctypes.c_double]

    return mp_fit


class MPFit(daoFitC.MultiFitterBase):
    """
    Base class for multi-plane fitting. While technically a sub-class of 
    daoFitC.MultiFitterBase this overrides essentially all of the base 
    class methods.
    """
    def __init__(self, independent_heights = None, psf_objects = None, scmos_cals = None, **kwds):
        super(MPFit, self).__init__(**kwds)

        self.clib = loadMPFitC()
        self.independent_heights = independent_heights
        self.n_channels = len(psf_objects)
        self.psf_objects = psf_objects
        self.scmos_cals = []

    def cleanup(self, spacing = "  ", verbose = True):
        if self.mfit is not None:
            if verbose:
                for i in range(self.n_channels):
                    print("Channel", i)
                    daoFitC.printFittingInfo(self.mfit.contents.fit_data[i],
                                             spacing = spacing)
            print()
            print(spacing, "{0:0d} fitting iterations.".format(self.iterations))
            print(spacing, "{0:.1f} fitting iterations/channel.".format(float(self.iterations)/float(self.n_channels)))
            self.clib.mpCleanup(self.mfit)

    def getFitImage(self, channel):
        """
        Get a fit image, i.e. f(x), an image created from drawing all of
        the current fits into a 2D array for the specified channel.
        """
        fit_image = numpy.zeros(self.im_shape)
        self.clib.mpGetFitImage(self.mfit, fit_image, channel)
        return fit_image

    def getIterations(self):
        """
        Update iterations and reset C library counter. The idea anyway
        is that the Python counter won't overflow where as the C counter
        might, particularly on a 32 bit system.
        """
        for i in range(self.n_channels):
            fit_data = self.mfit.contents.fit_data[i]
            self.iterations += fit_data.contents.n_iterations
            fit_data.contents.n_iterations = 0
        return self.iterations

    def getGoodPeaks(self, peaks, threshold):
        #
        # FIXME: We'd like to have a minimum height threshold, but the
        #        threshold parameter is a relative value not an absolute.
        #        For now we are just rejecting ERROR peaks.
        #
        if (peaks.size > 0):
            #
            # In the C library, the peaks that represent a single object
            # in multiple channels will all have the same status index
            # and height, so we can filter here with some confidence that
            # we'll keep or eliminate all the peaks in a particular group
            # at the same time. We attempt to enforce this with an assertion
            # statement, but it is an imperfect test as we could eliminate
            # 1 peak from each of 4 groups and still pass, even though this
            # would completely mess up the analysis.
            #
            status_index = utilC.getStatusIndex()
            height_index = utilC.getHeightIndex()

            #mask = (peaks[:,status_index] != 2.0) & (peaks[:,height_index] > min_height)
            mask = (peaks[:,status_index] != 2.0)
            if self.verbose:
                print(" ", numpy.sum(mask), "were good out of", peaks.shape[0])

            #
            # For debugging.
            #
            if False:
                xw_index = utilC.getXWidthIndex()
                yw_index = utilC.getYWidthIndex()
                for i in range(peaks.shape[0]):
                    if (peaks[i,xw_index] != peaks[i,yw_index]):
                        print("bad peak width detected", i, peaks[i,xw_index], peaks[i,yw_index])

            #
            # Debugging check that the peak status markings are actually
            # in sync.
            #
            if True:
                n_peaks = int(peaks.shape[0] / self.n_channels)
                not_bad = True
                for i in range(n_peaks):
                    if (mask[i] != mask[i+n_peaks]):
                        print("Problem detected with peak", i)
                        print("  ", peaks[i,:])
                        print("  ", peaks[i+n_peaks,:])
                        not_bad = False
                assert not_bad

            masked_peaks = peaks[mask,:]

            # The number of peaks should always be a multiple of the number of channels.
            assert((masked_peaks.shape[0] % self.n_channels) == 0)
                
            return masked_peaks

        else:
            return peaks

    def getResults(self, input_peaks_shape):
        fit_peaks = numpy.ascontiguousarray(numpy.zeros(input_peaks_shape))
        self.clib.mpGetResults(self.mfit, fit_peaks)
        return fit_peaks

    def getUnconverged(self):
        return self.clib.mpGetUnconverged(self.mfit)

    def initializeC(self, variance):
        self.im_shape = variance.shape

        # Get fitting structure.
        self.mfit = self.clib.mpInitialize(numpy.ascontiguousarray(self.clamp),
                                           self.default_tol,
                                           self.n_channels,
                                           self.independent_heights,
                                           variance.shape[1],
                                           variance.shape[0])
    
    def iterate(self):
        self.clib.mpIterateLM(self.mfit)
        #self.clib.mpIterateOriginal(self.mfit)

    def newImage(self, image, channel):
        if (image.shape[0] != self.im_shape[0]) or (image.shape[1] != self.im_shape[1]):
            raise daoFitC.MultiFitterException("Current image shape and the original image shape are not the same.")

        self.clib.mpNewImage(self.mfit, image, channel)

    def newPeaks(self, peaks):
        assert ((peaks.shape[0]%self.n_channels) == 0)

        n_peaks = int(peaks.shape[0]/self.n_channels)
        self.clib.mpNewPeaks(self.mfit,
                             numpy.ascontiguousarray(peaks),
                             n_peaks)

    def rescaleZ(self, peaks):
        return peaks

    def setMapping(self, xt_0toN, yt_0toN, xt_Nto0, yt_Nto0):
        self.clib.mpSetTransforms(self.mfit,
                                  numpy.ascontiguousarray(xt_0toN),
                                  numpy.ascontiguousarray(yt_0toN),
                                  numpy.ascontiguousarray(xt_Nto0),
                                  numpy.ascontiguousarray(yt_Nto0))

    def setVariance(self, variance, channel):
        #
        # This is a little different than 3D-DAOSTORM, sCMOS and Spliner
        # because this is the first thing that will be called with an array
        # that is expected to have the same size of the images. So we initialize
        # the C library now, rather than in newImage().
        #
        if self.mfit is None:
            self.initializeC(variance)
        else:
            if (variance.shape[0] != self.im_shape[0]) or (variance.shape[1] != self.im_shape[1]):
                raise daoFitC.MultiFitterException("Current variance shape and the original variance shape are not the same.")

    def setWeights(self, weights, z_offset, z_scale):
        """
        Initialize weights. These are used to weight the per channel parameter
        update values based on the localizations z value. The idea is that
        at any particular z value some channels will contribute more information
        than others to the value of each parameter.

        The C library expects these arrays to be indexed by z value, then channel.
        """
        assert(weights["h"].shape[1] == self.n_channels), "Incorrect shape for weights array."

        w_bg = numpy.ascontiguousarray(weights["bg"])
        w_h = numpy.ascontiguousarray(weights["h"])
        w_x = numpy.ascontiguousarray(weights["x"])
        w_y = numpy.ascontiguousarray(weights["y"])
        w_z = numpy.ascontiguousarray(weights["z"])

        # Check for negative or zero values in the weights.
        for w in [w_bg, w_h, w_x, w_y, w_z]:
            mask = (w > 0.0)
            assert(numpy.count_nonzero(mask) == w.size)

        self.clib.mpSetWeights(self.mfit, w_bg, w_h, w_x, w_y, w_z, w_bg.shape[0])
        self.clib.mpSetWeightsIndexing(self.mfit, z_offset, z_scale)


class MPPupilFnFit(MPFit):
    """
    The basic idea is that we are going to use the functionality from PupilFn
    to do most of the work. We will have one pupilFit C structure per image
    plane / channel.
    """
    def __init__(self, **kwds):
        super(MPPupilFnFit, self).__init__(**kwds)

        self.clamp = numpy.array([1.0,  # Height (Note: This is relative to the initial guess).
                                  1.0,  # x position
                                  0.3,  # width in x (Note: Not relevant for this fitter).
                                  1.0,  # y position
                                  0.3,  # width in y (Note: Not relevant for this fitter).
                                  1.0,  # background (Note: This is relative to the initial guess).
                                  1.0]) # z position

    def setVariance(self, variance, channel):
        super(MPSplineFit, self).setVariance(variance, channel)
        
        # This where the differentation in which type of fitter to use happens.
        zmax = self.psf_objects[0].getZMax() * 1.0e-3
        zmin = self.psf_objects[0].getZMin() * 1.0e-3
        self.clib.mpInitializePupilFnChannel(self.mfit,
                                             self.psf_objects[channel].getCPointer(),
                                             variance,
                                             zmin,
                                             zmax,
                                             channel)

    def setWeights(self, weights):
        if weights is None:
            weights = {"bg" : numpy.ones((1, self.n_channels))/float(self.n_channels),
                       "h" : numpy.ones((1, self.n_channels))/float(self.n_channels),
                       "x" : numpy.ones((1, self.n_channels))/float(self.n_channels),
                       "y" : numpy.ones((1, self.n_channels))/float(self.n_channels),
                       "z" : numpy.ones((1, self.n_channels))/float(self.n_channels)}
            super(MPPupilFnFit, self).setWeights(weights, 0.0, 0.0)

        else:
            z_offset = self.psf_objects[0].getZMin() * 1.0e-3
            z_scale = float(weights["bg"].shape[0])/((self.psf_objects[0].getZMax() - self.psf_objects[0].getZMin() + 1.0e-9) * 1.0e-3)
            super(MPPupilFnFit, self).setWeights(weights, z_offset, z_scale)
            
    
class MPSplineFit(MPFit):
    """
    The basic idea is that we are going to use the functionality from Spliner
    to do most of the work. We will have one splineFit C structure per image
    plane / channel.
    """
    def __init__(self, **kwds):
        super(MPSplineFit, self).__init__(**kwds)

        # Clamp parameters.
        #
        if True:
            self.clamp = numpy.array([1.0,  # Height (Note: This is relative to the initial guess).
                                      1.0,  # x position
                                      0.3,  # width in x (Note: Not relevant for this fitter).
                                      1.0,  # y position
                                      0.3,  # width in y (Note: Not relevant for this fitter).
                                      1.0,  # background (Note: This is relative to the initial guess).
                                      0.5 * self.psf_objects[0].getSplineSize()]) # z position (in spline size units).
        else:
            self.clamp = numpy.array([1.0,  # Height (Note: This is relative to the initial guess).
                                      1.0,  # x position
                                      0.3,  # width in x (Note: Not relevant for this fitter).
                                      1.0,  # y position
                                      0.3,  # width in y (Note: Not relevant for this fitter).
                                      1.0,  # background (Note: This is relative to the initial guess).
                                      1.0]) # z position (in spline size units).

    def rescaleZ(self, peaks):
        z_index = utilC.getZCenterIndex()
        peaks[:,z_index] = self.psf_objects[0].rescaleZ(peaks[:,z_index])
        return peaks
    
    def setVariance(self, variance, channel):
        super(MPSplineFit, self).setVariance(variance, channel)
        
        # This where the differentation in which type of fitter to use happens.
        self.clib.mpInitializeSplineChannel(self.mfit,
                                            self.psf_objects[channel].getCPointer(),
                                            variance,
                                            channel)

    def setWeights(self, weights):
        if weights is None:
            weights = {"bg" : numpy.ones((1, self.n_channels))/float(self.n_channels),
                       "h" : numpy.ones((1, self.n_channels))/float(self.n_channels),
                       "x" : numpy.ones((1, self.n_channels))/float(self.n_channels),
                       "y" : numpy.ones((1, self.n_channels))/float(self.n_channels),
                       "z" : numpy.ones((1, self.n_channels))/float(self.n_channels)}
            super(MPSplineFit, self).setWeights(weights, 0.0, 0.0)

        else:
            assert (weights["bg"].shape[0] == self.psf_objects[0].getSplineSize()), "Incorrect shape for weights array."
            super(MPSplineFit, self).setWeights(weights, 0.0, 1.0)
            
