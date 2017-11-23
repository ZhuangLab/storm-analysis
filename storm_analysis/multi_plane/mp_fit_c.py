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

                ('independent_heights', ctypes.c_int),

                ('n_channels', ctypes.c_int),
                ('n_weights', ctypes.c_int),
                
                ('w_z_offset', ctypes.c_double),
                ('w_z_scale', ctypes.c_double),
                
                ('tolerance', ctypes.c_double),

                ('zmin', ctypes.c_double),
                ('zmax', ctypes.c_double),
                
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
                
                ('fn_cleanup', ctypes.c_void_p),
                ('fn_newpeaks', ctypes.c_void_p),
                ('fn_peak_xi_yi', ctypes.c_void_p),
                ('fn_update', ctypes.c_void_p),
                ('fn_z_range', ctypes.c_void_p)]


def loadMPFitC():
    mp_fit = loadclib.loadCLibrary("storm_analysis.multi_plane", "mp_fit")
    
    # These are from sa_library/multi_fit.c
    mp_fit.mFitGetFitImage.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]
    
    mp_fit.mFitGetNError.argtypes = [ctypes.c_void_p]
    mp_fit.mFitGetNError.restype = ctypes.c_int
    
    mp_fit.mFitGetPeakPropertyDouble.argtypes = [ctypes.c_void_p,
                                                 ndpointer(dtype=numpy.float64),
                                                 ctypes.c_char_p]

    mp_fit.mFitGetPeakPropertyInt.argtypes = [ctypes.c_void_p,
                                              ndpointer(dtype=numpy.int32),
                                              ctypes.c_char_p]
    
    mp_fit.mFitGetResidual.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=numpy.float64)]
    
    mp_fit.mFitGetUnconverged.argtypes = [ctypes.c_void_p]
    mp_fit.mFitGetUnconverged.restype = ctypes.c_int

    mp_fit.mFitIterateLM.argtypes = [ctypes.c_void_p]
    mp_fit.mFitIterateOriginal.argtypes = [ctypes.c_void_p]

    mp_fit.mFitNewBackground.argtypes = [ctypes.c_void_p,
                                         ndpointer(dtype=numpy.float64)]
    
    mp_fit.mFitNewImage.argtypes = [ctypes.c_void_p,
                                    ndpointer(dtype=numpy.float64)]

    mp_fit.mFitRemoveErrorPeaks.argtypes = [ctypes.c_void_p]

    mp_fit.mFitRemoveRunningPeaks.argtypes = [ctypes.c_void_p]

    mp_fit.mFitSetPeakStatus.argtypes = [ctypes.c_void_p,
                                         ndpointer(dtype=numpy.int32)]
    
    # From multi_plane/mp_fit.c
    mp_fit.mpCleanup.argtypes = [ctypes.c_void_p]
    
    mp_fit.mpInitialize.argtypes = [ndpointer(dtype=numpy.float64),
                                    ctypes.c_double,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]
    mp_fit.mpInitialize.restype = ctypes.POINTER(mpFitData)

    mp_fit.mpInitializePSFFFTChannel.argtypes = [ctypes.c_void_p,
                                                 ctypes.c_void_p,
                                                 ndpointer(dtype=numpy.float64),
                                                 ctypes.c_int]
    
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

    mp_fit.mpNewPeaks.argtypes = [ctypes.c_void_p,
                                  ndpointer(dtype=numpy.float64),
                                  ctypes.c_char_p,
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


class MPFit(daoFitC.MultiFitterArbitraryPSF):
    """
    Base class for multi-plane fitting.
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
            print(spacing, self.n_proximity, "peaks lost to proximity.")
            print(spacing, "{0:0d} fitting iterations.".format(self.iterations))
            print(spacing, "{0:.1f} fitting iterations/channel.".format(float(self.iterations)/float(self.n_channels)))
            self.clib.mpCleanup(self.mfit)

    def getFitImage(self):
        """
        Return the fit images for each channel.
        """
        fit_images = []
        for i in range(self.n_channels):
            fit_image = numpy.zeros(self.im_shape)
            self.clib.mFitGetFitImage(self.mfit.contents.fit_data[i], fit_image)
            fit_images.append(fit_image)
        return fit_images

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

    def getNError(self, channel = 0):
        """
        Return the number of peaks in the error state.
        """
        return self.clib.mFitGetNError(self.mfit.contents.fit_data[channel])

    def getNFit(self, channel = 0):
        """
        Return the current number of peaks.
        """
        return self.mfit.contents.fit_data[channel].contents.nfit

    def getNFitMax(self, channel = 0):
        """
        Return the current number of peaks.
        """
        return self.mfit.contents.fit_data[channel].contents.max_nfit

    def getPeakProperty(self, p_name, channel = 0):
        """
        Return a numpy array containing the requested property.
        """
        if not p_name in self.peak_properties:
            raise MultiFitterException("No such property '" + p_name + "'")

        if(self.peak_properties[p_name] == "float"):
            values = numpy.ascontiguousarray(numpy.zeros(self.getNFit(), dtype = numpy.float64))
            self.clib.mFitGetPeakPropertyDouble(self.mfit.contents.fit_data[channel],
                                                values,
                                                ctypes.c_char_p(p_name.encode()))
            return values
        
        elif(self.peak_properties[p_name] == "int"):
            values = numpy.ascontiguousarray(numpy.zeros(self.getNFit(), dtype = numpy.int32))
            self.clib.mFitGetPeakPropertyInt(self.mfit.contents.fit_data[channel],
                                             values,
                                             ctypes.c_char_p(p_name.encode()))
            return values

    def getSize(self):
        return self.psf_objects[0].getSize()
        
    def getUnconverged(self, channel = 0):
        return self.clib.mFitGetUnconverged(self.mfit.contents.fit_data[channel])

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

        # FIXME: This option no longer works due to bit rot in the C library.
        #self.clib.mpIterateOriginal(self.mfit)

    def newBackground(self, background):
        """
        background - a list of background estimates of length n_channels.
        """
        if (len(background) != self.n_channels):
            raise daoFitC.MultiFitterException("Number of images does not match number of channels.")
        
        for i in range(self.n_channels):
            if (background[i].shape[0] != self.im_shape[0]) or (background[i].shape[1] != self.im_shape[1]):
                raise daoFitC.MultiFitterException("Background image shape and the original image shape are not the same.")

            self.clib.mFitNewBackground(self.mfit.contents.fit_data[i],
                                        numpy.ascontiguousarray(background[i], dtype = numpy.float64))

    def newImage(self, image):
        """
        image - a list of images of length n_channels.
        """
        if (len(image) != self.n_channels):
            raise daoFitC.MultiFitterException("Number of images does not match number of channels.")
        
        for i in range(self.n_channels):
            if (image[i].shape[0] != self.im_shape[0]) or (image[i].shape[1] != self.im_shape[1]):
                raise daoFitC.MultiFitterException("Current image shape and the original image shape are not the same.")

            self.clib.mFitNewImage(self.mfit.contents.fit_data[i],
                                   numpy.ascontiguousarray(image[i], dtype = numpy.float64))            

    def newPeaks(self, peaks, peaks_type):
        c_peaks = self.formatPeaks(peaks, peaks_type)
        self.clib.mpNewPeaks(self.mfit,
                             c_peaks,
                             ctypes.c_char_p(peaks_type.encode()),
                             c_peaks.shape[0])

    def removeErrorPeaks(self, check = True):
        for i in range(self.n_channels):
            self.clib.mFitRemoveErrorPeaks(self.mfit.contents.fit_data[i])

        if check:
            for i in range(1, self.n_channels):
                assert(self.getNFit(0) == self.getNFit(i))

    def removeRunningPeaks(self, check = True):
        for i in range(self.n_channels):
            self.clib.mFitRemoveRunningPeaks(self.mfit.contents.fit_data[i])

        if check:
            for i in range(1, self.n_channels):
                assert(self.getNFit(0) == self.getNFit(i))                

    def rescaleZ(self, z):
        return z

    def setMapping(self, xt_0toN, yt_0toN, xt_Nto0, yt_Nto0):
        self.clib.mpSetTransforms(self.mfit,
                                  numpy.ascontiguousarray(xt_0toN),
                                  numpy.ascontiguousarray(yt_0toN),
                                  numpy.ascontiguousarray(xt_Nto0),
                                  numpy.ascontiguousarray(yt_Nto0))

    def setPeakStatus(self, status):
        """
        Set the status (RUNNING, CONVERGED, ERROR) of the peaks in the C library.
        """
        assert (status.size == self.getNFit())
        for i in range(self.n_channels):
            self.clib.mFitSetPeakStatus(self.mfit.contents.fit_data[i],
                                        numpy.ascontiguousarray(status, dtype = numpy.int32))
        
    def setVariance(self, variance, channel):
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
        print("weights z scaling - offset: {0:.3f} scale: {1:.3f}".format(z_offset, z_scale))
        self.clib.mpSetWeightsIndexing(self.mfit, z_offset, z_scale)


class MPPSFFnFit(MPFit):
    """
    The basic idea is that we are going to use the functionality from PSF FFT
    to do most of the work. We will have one psfFFTFit C structure per image
    plane / channel.
    """
    def __init__(self, **kwds):
        super(MPPSFFnFit, self).__init__(**kwds)

        self.clamp = numpy.array([1.0,  # Height (Note: This is relative to the initial guess).
                                  1.0,  # x position
                                  0.3,  # width in x (Note: Not relevant for this fitter).
                                  1.0,  # y position
                                  0.3,  # width in y (Note: Not relevant for this fitter).
                                  1.0,  # background (Note: This is relative to the initial guess).
                                  0.5 * self.psf_objects[0].getZSize()]) # z position (in FFT size units).

    def rescaleZ(self, z):
        return self.psf_objects[0].rescaleZ(z)
    
    def setVariance(self, variance, channel):
        super(MPPSFFnFit, self).setVariance(variance, channel)
        
        # This where the differentation in which type of fitter to use happens.
        zmax = self.psf_objects[0].getZMax() * 1.0e-3
        zmin = self.psf_objects[0].getZMin() * 1.0e-3
        self.clib.mpInitializePSFFFTChannel(self.mfit,
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
            super(MPSFFnFit, self).setWeights(weights, 0.0, 0.0)

        else:
            zmax = self.psf_objects[0].getZMax() * 1.0e-3
            zmin = self.psf_objects[0].getZMin() * 1.0e-3
            z_offset = -0.5*float(self.psf_objects[0].getZSize())
            z_scale = float(weights["bg"].shape[0])/float(self.psf_objects[0].getZSize())
            super(MPPSFFnFit, self).setWeights(weights, z_offset, z_scale)

            
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
        super(MPPupilFnFit, self).setVariance(variance, channel)
        
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
            zmax = self.psf_objects[0].getZMax() * 1.0e-3
            zmin = self.psf_objects[0].getZMin() * 1.0e-3
            z_offset = zmin
            z_scale = float(weights["bg"].shape[0])/(zmax - zmin + 1.0e-12)
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

    def rescaleZ(self, z):
        return self.psf_objects[0].rescaleZ(z)
    
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
            assert (weights["bg"].shape[0] == (self.psf_objects[0].getSplineSize() + 1)), "Incorrect shape for weights array."
            super(MPSplineFit, self).setWeights(weights, 0.0, 1.0)
            
