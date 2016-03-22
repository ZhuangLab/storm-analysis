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

    objects = numpy.zeros((1, 5))

    # PSF of shape 1.
    objects[0,:] = [x_size/2, y_size/2, 1.0, 1.0, 1.0]
    psf1 = dg.drawGaussians((x_size, y_size), objects)
    psf1 = psf1/numpy.sum(psf1)
    flt1 = MatchedFilter(psf1)

    # PSF of shape 2.
    objects[0,:] = [x_size/2, y_size/2, 1.0, 2.0, 0.7]
    psf2 = dg.drawGaussians((x_size, y_size), objects)
    psf2 = psf2/numpy.sum(psf2)
    flt2 = MatchedFilter(psf2)

    result1 = flt1.convolve(10000.0 * psf2)
    result2 = flt2.convolve(10000.0 * psf2)

    print "Match to 1:", numpy.max(result1)
    print "Match to 2:", numpy.max(result2)

    with tifffile.TiffWriter("result.tif") as tif:
        tif.save(result1.astype(numpy.uint16))
        tif.save(result2.astype(numpy.uint16))

