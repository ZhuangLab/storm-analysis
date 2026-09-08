#!/usr/bin/env python
"""
Draws guassians onto a user supplied image.

Hazen 01/16
"""

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.loadclib as loadclib

drawgauss = loadclib.loadCLibrary("draw_gaussians")

drawgauss.drawGaussians.argtypes = [ndpointer(dtype = numpy.float64),
                                    ndpointer(dtype = numpy.float64),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]

def cDrawGaussians(image, objects, resolution):
    """
    Draws Gaussian's on image inplace at the requested resolution.
    
    image - A 2D C contiguous numpy.float64 array.
    objects - (N,5) numpy array, [[x1,y1,h1,sx1,sy1], [x2,y2,h2,sx2,sy2], ..]
    resolution - Number of sub-pixels in X/Y to sample over.
    """
    assert (image.dtype == numpy.float64), "Image type is not numpy.float64"
    assert (image.flags['C_CONTIGUOUS']), "Image is not C contiguous."
    assert (len(image.shape) == 2), "Image is not a 2D array."
    c_objects = numpy.ascontiguousarray(objects, dtype = numpy.float64)
    drawgauss.drawGaussians(image,
                            c_objects,
                            image.shape[1],
                            image.shape[0],
                            objects.shape[0],
                            resolution)

def drawGaussians(size, objects, background = 0.0, res = 1):
    """
    Creates an image, then draws Gaussians on it at the requested resolution.

    size - List specifying the image dimensions, [nx,ny].
    objects - (N,5) numpy array, [[x1,y1,h1,sx1,sy1], [x2,y2,h2,sx2,sy2], ..]
    background - A constant background term.
    res - Number of sub-pixels in X/Y to sample over.
    """
    assert (len(size)==2), "Image size must be (nx,ny)."
    image = numpy.ascontiguousarray(background * numpy.ones(size), dtype = numpy.float64)
    cDrawGaussians(image, objects, res)
    return image

def drawGaussiansXY(size, x, y, height = 1.0, sigma = 1.0, background = 0.0, res = 1):
    """
    Creates an image, then draws Gaussians at the requested x,y locations.

    size - List specifying the image dimensions, (nx,ny).
    x - Numpy array of Gaussian x positions.
    y - Numpy array of Gaussian y positions.
    height - Height to use for all of the Gaussians.
    sigma - Sigma to use for all of the Gaussians.
    background - A constant background term.
    res - Number of sub-pixels in X/Y to sample over.
    """
    assert (len(size)==2), "Image size must be (nx,ny)."
    image = background * numpy.ones(size)
    np = x.shape[0]
    h = numpy.ones(np) * float(height)
    sx = numpy.ones(np) * float(sigma)
    sy = numpy.ones(np) * float(sigma)
    objects = numpy.concatenate((x[:,None],
                                 y[:,None],
                                 h[:,None],
                                 sx[:,None],
                                 sy[:,None]),
                                axis = 1)
    cDrawGaussians(image, objects, res)
    return image

def drawGaussiansXYZ(size, x, y, z, height = 1.0, sigma = 1.0, background = 0.0, res = 1):
    """
    Creates an image, then draws Gaussians at the requested x,y,z locations.

    size - List specifying the image dimensions, (nz, nx, ny).
    x - Numpy array of Gaussian x positions.
    y - Numpy array of Gaussian y positions.
    z - Numpy array of Gaussian z positions.
    height - Height to use for all of the Gaussians.
    sigma - Sigma to use for all of the Gaussians.
    background - A constant background term.
    res - Number of sub-pixels in X/Y to sample over.
    """
    assert (len(size)==3), "Image size must be (nz,nx,ny)."
    image = background * numpy.ones(size)
    np = x.shape[0]
    sx = numpy.ones(np) * float(sigma)
    sy = numpy.ones(np) * float(sigma)
    for i in range(size[0]):
        dz = i - z
        h = numpy.ones(np) * float(height) * numpy.exp(-dz*dz/(2.0*sigma*sigma))
        objects = numpy.concatenate((x[:,None],
                                     y[:,None],
                                     h[:,None],
                                     sx[:,None],
                                     sy[:,None]),
                                    axis = 1)
        cDrawGaussians(image[i,:,:], objects, res)
    return image

def drawGaussiansXYOnImage(image, x, y, height = 1.0, sigma = 1.0, res = 1):
    """
    Draws Gaussians at the requested x,y locations on the image inplace.

    size - List specifying the image dimensions, (nx,ny).
    x - Numpy array of Gaussian x positions.
    y - Numpy array of Gaussian y positions.
    height - Height to use for all of the Gaussians.
    sigma - Sigma to use for all of the Gaussians.
    background - A constant background term.
    res - Number of sub-pixels in X/Y to sample over.
    """
    np = x.shape[0]
    h = numpy.ones(np) * float(height)
    sx = numpy.ones(np) * float(sigma)
    sy = numpy.ones(np) * float(sigma)
    objects = numpy.concatenate((x[:,None],
                                 y[:,None],
                                 h[:,None],
                                 sx[:,None],
                                 sy[:,None]),
                                axis = 1)
    cDrawGaussians(image, objects, res)
    return image


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
