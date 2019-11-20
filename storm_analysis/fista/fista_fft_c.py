#!/usr/bin/env python
"""
Simple Python interface to fista_fft.c

Hazen 2/16
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer

import storm_analysis.sa_library.cs_algorithm as csAlgorithm
import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.recenter_psf as recenterPSF

import storm_analysis.fista.fista_3d as fista3D

fista_fft = loadclib.loadCLibrary("fista_fft")

# C interface definition
fista_fft.cleanup.argtypes = [ctypes.c_void_p]

fista_fft.dwlsError.argtypes = [ctypes.c_void_p]
fista_fft.dwlsError.restype = ctypes.c_double

fista_fft.getAx.argtypes = [ctypes.c_void_p,
                            ndpointer(dtype=numpy.float64)]

fista_fft.getXVector.argtypes = [ctypes.c_void_p,
                                 ndpointer(dtype=numpy.float64)]

fista_fft.initialize2D.argtypes = [ndpointer(dtype=numpy.float64),
                                   ctypes.c_double,
                                   ctypes.c_int,
                                   ctypes.c_int,
                                   ctypes.c_int]
fista_fft.initialize2D.restype = ctypes.c_void_p

fista_fft.initialize3D.argtypes = [ndpointer(dtype=numpy.float64),
                                   ctypes.c_double,
                                   ctypes.c_int,
                                   ctypes.c_int,
                                   ctypes.c_int,
                                   ctypes.c_int]
fista_fft.initialize3D.restype = ctypes.c_void_p

fista_fft.iterate.argtypes = [ctypes.c_void_p,
                              ctypes.c_double]

fista_fft.l1Error.argtypes = [ctypes.c_void_p]
fista_fft.l1Error.restype = ctypes.c_double

fista_fft.l2Error.argtypes = [ctypes.c_void_p]
fista_fft.l2Error.restype = ctypes.c_double

fista_fft.newImage.argtypes = [ctypes.c_void_p,
                               ndpointer(dtype=numpy.float64),
                               ndpointer(dtype=numpy.float64)]

fista_fft.run.argtypes = [ctypes.c_void_p,
                          ctypes.c_double,
                          ctypes.c_int]


class FISTAFFTException(Exception):
    pass


class FISTA(csAlgorithm.CSAlgorithm):

    def __init__(self, psfs, timestep, dwls = False, **kwds):
        """
        Psfs is an array of psfs for different image planes (nx, ny, nz).
        
        The PSFs must be the same size as the image that will get analyzed.
        """
        super(FISTA, self).__init__(**kwds)
                
        self.dwls = dwls
        self.shape = psfs.shape

        c_psfs = numpy.zeros(self.shape)
        for i in range(self.shape[2]):
            c_psfs[:,:,i] = recenterPSF.recenterPSF(psfs[:,:,i])
        c_psfs = numpy.ascontiguousarray(c_psfs, dtype = numpy.float)
        self.c_fista = fista_fft.initialize3D(c_psfs, timestep, self.shape[0], self.shape[1], self.shape[2], self.dwls)

    def cleanup(self):
        fista_fft.cleanup(self.c_fista)
        self.c_fista = None

    def dwlsError(self):
        return fista_fft.dwlsError(self.c_fista)

    def getAx(self):
        ax = numpy.ascontiguousarray(numpy.zeros((self.shape[0], self.shape[1]), dtype = numpy.float))
        fista_fft.getAx(self.c_fista, ax)
        return ax
        
    def getXVector(self):
        x_vector = numpy.ascontiguousarray(numpy.zeros(self.shape, dtype = numpy.float))
        fista_fft.getXVector(self.c_fista, x_vector)
        return x_vector
    
    def iterate(self, f_lambda):
        fista_fft.iterate(self.c_fista, f_lambda)

    def l1Error(self):
        return fista_fft.l1Error(self.c_fista)

    def l2Error(self):
        return fista_fft.l2Error(self.c_fista)
    
    def newImage(self, image, background):
        """
        image - The image to deconvolve, includes the background estimate.
        background - An estimate of the background in the image.
        """
        if (image.shape[0] != self.shape[0]) or (image.shape[1] != self.shape[1]):
            tmp_str = " ".join("Image and PSF are not the same size", str(image.shape), str(self.shape[:2]))
            raise FISTAFFTException(tmp_str)

        if (background.shape[0] != self.shape[0]) or (background.shape[1] != self.shape[1]):
            tmp_str = " ".join("Background and PSF are not the same size", str(background.shape), str(self.shape[:2]))
            raise FISTAFFTException(tmp_str)

        c_image = numpy.ascontiguousarray(image, dtype = numpy.float)
        c_background = numpy.ascontiguousarray(background, dtype = numpy.float)
        fista_fft.newImage(self.c_fista, c_image, c_background)

    def run(self, f_lamba, iterations):
        fista_fft.run(self.c_fista, f_lambda, iterations)


#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
