#!/usr/bin/env python
"""
Simple Python interface to psf_fft.c

Hazen 2/17
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

import storm_analysis.sa_library.loadclib as loadclib

psf_fft = loadclib.loadCLibrary("storm_analysis.fft_fitting", "psf_fft")

psf_fft.cleanup.argtypes = [ctypes.c_void_p]

psf_fft.getPSF.argtypes = [ctypes.c_void_p,
                           ndpointer(dtype = numpy.float64),
                           ctypes.c_double,
                           ctypes.c_double,
                           ctypes.c_double]

psf_fft.initialize.argtypes = [ndpointer(dtype = numpy.float64),
                               ctypes.c_int,
                               ctypes.c_int,
                               ctypes.c_int]
psf_fft.initialize.restype = ctypes.c_void_p


class PSFFFTException(Exception):

    def __init__(self, message):
        Exception.__init__(self, message)

        
class PSFFFT(object):

    def __init__(self, psf):
        self.psf_shape = psf.shape

        c_psf = numpy.ascontiguousarray(psf, dtype = numpy.float64)
        self.pfft = psf_fft.initialize(c_psf,
                                       self.psf_shape[0],
                                       self.psf_shape[1],
                                       self.psf_shape[2])

    def cleanup(self):
        psf_fft.cleanup(self.pfft)

    def getPSF(self, dx, dy, dz):
        psf = numpy.zeros((self.psf_shape[1], self.psf_shape[2]), dtype = numpy.float64)
        psf_fft.getPSF(self.pfft, psf, dx, dy, dz)
        return psf


if (__name__ == "__main__"):

    psf = numpy.array([[[0, 0, 0],
                        [0, 1, 0],
                        [0, 0, 0]],
                       [[0, 1, 0],
                        [1, 2, 1],
                        [0, 1, 0]],
                       [[0, 0, 0],
                        [0, 1, 0],
                        [0, 0, 0]]])

    psf = numpy.array([[[0, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 0]],
                       [[0, 1, 0, 0],
                        [1, 2, 1, 0],
                        [0, 1, 0, 0]],
                       [[0, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 0]]])
    
    pfft = PSFFFT(psf)

    print(pfft.getPSF(0,0,1))
    
    pfft.cleanup()

#
# The MIT License
#
# Copyright (c) 2017 Zhuang Lab, Harvard University
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
