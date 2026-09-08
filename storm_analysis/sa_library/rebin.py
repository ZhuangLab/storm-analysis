#!/usr/bin/env python
"""
Rebinning of arrays.

Hazen 2/15
"""

import numpy
import scipy
import scipy.fftpack

#
# This comes directly from here:
#
# http://scipy-cookbook.readthedocs.io/items/Rebinning.html
#
# Copyright 2015, Various authors. Revision 28f8936a.
# 
def downSample(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = numpy.asarray(shape)/numpy.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],int(factor[%d]),'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    return eval(''.join(evList))

def upSampleFFT(image, factor):
    """
    Upsample using a FFT (high frequencies are set to zero).
    """
    xsize = image.shape[0]*factor
    ysize = image.shape[1]*factor

    new_image_fft = numpy.zeros((xsize, ysize),dtype=numpy.complex)
    image_fft = scipy.fftpack.fft2(image)
    half_x = int(image.shape[0]/2)
    half_y = int(image.shape[1]/2)
    new_image_fft[:half_x,:half_y] = image_fft[:half_x,:half_y]
    new_image_fft[-half_x:,:half_y] = image_fft[half_x:,:half_y]
    new_image_fft[:half_x,-half_y:] = image_fft[:half_x,half_y:]
    new_image_fft[-half_x:,-half_y:] = image_fft[half_x:,half_y:]

    new_image = numpy.real(scipy.fftpack.ifft2(new_image_fft))
    new_image[(new_image<0.0)] = 0.0
    return new_image

#
# The MIT License
#
# Copyright (c) 2015 Zhuang Lab, Harvard University
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
