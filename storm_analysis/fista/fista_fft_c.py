#!/usr/bin/python
#
# Simple Python interface to fista_fft.c
#
# Hazen 2/16
#

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import sa_library.loadclib as loadclib

import fista_3d as fista3D

fista_fft = loadclib.loadCLibrary(os.path.dirname(__file__), "fista_fft")

# C interface definition
fista_fft.getXVector.argtypes = [ndpointer(dtype=numpy.float64)]
fista_fft.initialize2D.argtypes = [ndpointer(dtype=numpy.float64),
                                   ctypes.c_double,
                                   ctypes.c_int]
fista_fft.initialize3D.argtypes = [ndpointer(dtype=numpy.float64),
                                   ctypes.c_double,
                                   ctypes.c_int,
                                   ctypes.c_int]
fista_fft.iterate.argtypes = [ctypes.c_double]
fista_fft.l1Error.restype = ctypes.c_double
fista_fft.l2Error.restype = ctypes.c_double
fista_fft.newImage.argtypes = [ndpointer(dtype=numpy.float64)]
fista_fft.run.argtypes = [ctypes.c_double,
                          ctypes.c_int]


#
# Note that while this might appear like a nice object, there are a lot of
# static variables in the fista_fft C library so in reality you should only
# be using one of these per process otherwise chaos will likely ensue.
#
class FISTA(object):

    #
    # For 2D psfs is an array of size (nx, nx).
    #
    # For 3D psfs is an array of psfs for different image planes (nx, nx, nz).
    #
    # The PSFs must be the same size as the image that will get analyzed.
    #
    def __init__(self, psfs, timestep):

        if (psfs.shape[0] != psfs.shape[1]):
            print "The PSF must be square (in X-Y)!"
            exit()
            
        self.shape = psfs.shape

        if (len(self.shape) == 2):
            c_psfs = numpy.ascontiguousarray(fista3D.recenterPSF(psfs), dtype = numpy.float)
            fista_fft.initialize2D(c_psfs, timestep, self.shape[0])
        else:
            c_psfs = numpy.zeros(self.shape)
            for i in range(self.shape[2]):
                c_psfs[:,:,i] = fista3D.recenterPSF(psfs[:,:,i])
            c_psfs = numpy.ascontiguousarray(c_psfs, dtype = numpy.float)
            fista_fft.initialize3D(c_psfs, timestep, self.shape[0], self.shape[2])

    def getXVector(self):
        x_vector = numpy.ascontiguousarray(numpy.zeros(self.shape, dtype = numpy.float))
        fista_fft.getXVector(x_vector)
        return x_vector
    
    def iterate(self, f_lambda):
        fista_fft.iterate(f_lambda)

    def l1Error(self):
        return fista_fft.l1Error()

    def l2Error(self):
        return fista_fft.l2Error()
    
    def newImage(self, image):
        if (image.shape[0] != self.shape[0]) or (image.shape[1] != self.shape[1]):
            print "Image and PSF are not the same size", image.shape, self.shape[:2]
            return
        c_image = numpy.ascontiguousarray(image, dtype = numpy.float)
        fista_fft.newImage(c_image)

    def run(self, f_lamba, iterations):
        fista_fft.run(f_lambda, iterations)        


if (__name__ == "__main__"):

    # Clever test here..
    pass


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
