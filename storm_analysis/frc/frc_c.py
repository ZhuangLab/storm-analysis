#!/usr/bin/python
#
# Simple Python interface to frc.c.
#
# Hazen 10/14
#

from ctypes import *
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

dll_path = os.path.dirname(__file__)
if not (dll_path == ""):
    dll_path += "/"

if(os.path.exists(dll_path + "frc.so")):
    frc_lib = cdll.LoadLibrary(dll_path + "frc.so")
else:
    frc_lib = cdll.LoadLibrary(dll_path + "frc.dll")

# Function specifications.
frc_lib.calc_frc.argtypes = [ndpointer(dtype=numpy.complex128),
                             ndpointer(dtype=numpy.complex128),
                             ndpointer(dtype=numpy.float64),                    
                             ndpointer(dtype=numpy.int32),
                             c_int,
                             c_int]

# Calculate FRC (this assumes that the images are square.
def frc(fft1, fft2):
    y_size = fft1.shape[0]
    x_size = fft1.shape[1]

    c_fft1 = numpy.ascontiguousarray(fft1).astype(numpy.complex128)
    c_fft2 = numpy.ascontiguousarray(fft2).astype(numpy.complex128)
    c_frc = numpy.zeros(x_size, dtype=numpy.float64)
    c_frc_counts = numpy.zeros(x_size, dtype=numpy.int32)
    frc_lib.calc_frc(c_fft1,
                     c_fft2,
                     c_frc,
                     c_frc_counts,
                     y_size,
                     x_size)

    max_q = min([x_size,y_size])/2
    return c_frc[:max_q], c_frc_counts[:max_q]


if (__name__ == "__main__"):
    test = numpy.ones((4,4))
    #test[1,1] = 0.0

    fft1 = numpy.fft.fftshift(numpy.fft.fft2(test))
    fft2 = numpy.fft.fftshift(numpy.fft.fft2(test))

    print fft1

    [results, counts] = frc(fft1, fft2)

    print results
    print counts

#
# The MIT License
#
# Copyright (c) 2014 Zhuang Lab, Harvard University
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
