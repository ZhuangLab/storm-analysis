#!/usr/bin/env python
#
# Simply Python interface to matched_filter.c
#
# Hazen 3/16
#

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

import sa_library.loadclib as loadclib
import sa_library.recenter_psf as recenterPSF

m_filter = loadclib.loadCLibrary(os.path.dirname(__file__), "matched_filter")

m_filter.cleanup.argtypes = [ctypes.c_void_p]
m_filter.convolve.argtypes = [ctypes.c_void_p,
                              ndpointer(dtype = numpy.float64),
                              ndpointer(dtype = numpy.float64)]
m_filter.initialize.argtypes = [ndpointer(dtype = numpy.float64),
                                ctypes.c_int,
                                ctypes.c_int]
m_filter.initialize.restype = ctypes.c_void_p


class MatchedFilterException(Exception):

    def __init__(self, message):
        Exception.__init__(self, message)

        
class MatchedFilter(object):

    def __init__(self, psf):
        self.psf_shape = psf.shape

        rc_psf = recenterPSF.recenterPSF(psf)

        self.mfilter = m_filter.initialize(rc_psf, rc_psf.shape[0], rc_psf.shape[1])

    def cleanup(self):
        m_filter.cleanup(self.mfilter)

    def convolve(self, image):
        if (image.shape[0] != self.psf_shape[0]) or (image.shape[1] != self.psf_shape[1]):
            raise MatchedFilterException("Image shape must match psf shape! " + str(image.shape) + " != " + str(self.psf_shape))
        
        image = numpy.ascontiguousarray(image, dtype = numpy.float64)
        result = numpy.zeros(self.psf_shape, dtype = numpy.float64)
        m_filter.convolve(self.mfilter, image, result)

        return result
        

if (__name__ == "__main__"):
    
    import tifffile

    import simulator.drawgaussians as dg

    x_size = 200
    y_size = 300

    # Create PSF and normalize area to 1.
    psf = dg.drawGaussiansXY((x_size, y_size),
                             numpy.array([x_size/2]),
                             numpy.array([y_size/2]))
    psf = psf/numpy.sum(psf)

    psf_16 = (1000.0 * psf).astype(numpy.uint16)
    tifffile.imsave("psf.tif", psf_16)

    flt = MatchedFilter(psf)
    
    # Make a fake image and convolve with the PSF.
    image = numpy.zeros((x_size, y_size))
    image[50,100] = 1000.0
    image[70,100] = 1000.0

    result = flt.convolve(image)
    print numpy.max(result)
    result = result.astype(numpy.uint16)
    tifffile.imsave("result.tif", result)

    flt.cleanup()

