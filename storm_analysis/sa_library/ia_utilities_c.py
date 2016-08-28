#!/usr/bin/python
#
# Simple Python interface to utilities.c.
#
# Hazen 6/11
#

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import sa_library.loadclib as loadclib

util = loadclib.loadCLibrary(os.path.dirname(__file__), "ia_utilities")

# C interface definition
util.findLocalMaxima.argtypes = [ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.int32),
                                 ndpointer(dtype=numpy.float64),
                                 ctypes.c_double,
                                 ctypes.c_double,
                                 ctypes.c_int,
                                 ctypes.c_int,
                                 ctypes.c_int,
                                 ctypes.c_int]
util.findLocalMaxima.restype = ctypes.c_int
util.getBackgroundIndex.restype = ctypes.c_int
util.getErrorIndex.restype = ctypes.c_int
util.getHeightIndex.restype = ctypes.c_int
util.getNPeakPar.restype = ctypes.c_int
util.getNResultsPar.restype = ctypes.c_int
util.getStatusIndex.restype = ctypes.c_int
util.getXCenterIndex.restype = ctypes.c_int
util.getXWidthIndex.restype = ctypes.c_int
util.getYCenterIndex.restype = ctypes.c_int
util.getYWidthIndex.restype = ctypes.c_int
util.getZCenterIndex.restype = ctypes.c_int
util.initializePeaks.argtypes = [ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.float64),
                                 ctypes.c_double,
                                 ctypes.c_double,
                                 ctypes.c_int,
                                 ctypes.c_int]
util.mergeNewPeaks.argtypes = [ndpointer(dtype=numpy.float64),
                               ndpointer(dtype=numpy.float64),
                               ndpointer(dtype=numpy.float64),
                               ctypes.c_double,
                               ctypes.c_double,
                               ctypes.c_int,
                               ctypes.c_int]
util.mergeNewPeaks.restype = ctypes.c_int
util.peakToPeakDist.argtypes = [ndpointer(dtype=numpy.float64),
                                ndpointer(dtype=numpy.float64),
                                ndpointer(dtype=numpy.float64),
                                ndpointer(dtype=numpy.float64),
                                ndpointer(dtype=numpy.float64),
                                ctypes.c_int, 
                                ctypes.c_int]
util.peakToPeakIndex.argtypes = [ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.int32),
                                 ctypes.c_int, 
                                 ctypes.c_int]
util.removeClosePeaks.argtypes = [ndpointer(dtype=numpy.float64),
                                  ndpointer(dtype=numpy.float64),
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_int]
util.removeClosePeaks.restype = ctypes.c_int
util.removeNeighbors.argtypes = [ndpointer(dtype=numpy.float64),
                                 ndpointer(dtype=numpy.float64),
                                 ctypes.c_double,
                                 ctypes.c_int]
util.removeNeighbors.restype = ctypes.c_int
util.smoothImage.argtypes = [ndpointer(dtype=numpy.float64),
                             ctypes.c_int]


# Return locations of local maxima
def findLocalMaxima(image, taken, threshold, radius, margin, maxpeaks = 10000):
    n_peak_par = getNResultsPar()
    image_c = numpy.ascontiguousarray(image)
    taken_c = numpy.ascontiguousarray(taken)
    peaks = numpy.ascontiguousarray(numpy.zeros(maxpeaks*n_peak_par))
    counts = util.findLocalMaxima(image_c,
                                  taken_c,
                                  peaks,
                                  threshold,
                                  radius,
                                  image.shape[1],
                                  image.shape[0],
                                  margin,
                                  maxpeaks)
    if (counts == maxpeaks):
        print "Found too many peaks..", counts
    peaks.resize(counts*n_peak_par)
    return [numpy.reshape(peaks, (-1,n_peak_par)), taken_c]

# Get the index of the background parameter.
def getBackgroundIndex():
    return util.getBackgroundIndex()

# Get the index of the error parameter.
def getErrorIndex():
    return util.getErrorIndex()

# Get the index of the height parameter.
def getHeightIndex():
    return util.getHeightIndex()

# Get the number of parameters in the peak fitting array.
def getNPeakPar():
    return util.getNPeakPar()

# Get the number of parameters in the results array.
def getNResultsPar():
    return util.getNResultsPar()

# Get the index of the status parameter
def getStatusIndex():
    return util.getStatusIndex()

# Get the index of the xcenter parameter
def getXCenterIndex():
    return util.getXCenterIndex()

# Get the index of the xwidth parameter
def getXWidthIndex():
    return util.getXWidthIndex()

# Get the index of the ycenter parameter
def getYCenterIndex():
    return util.getYCenterIndex()

# Get the index of the ywidth parameter
def getYWidthIndex():
    return util.getYWidthIndex()

# Get the index of the zcenter parameter
def getZCenterIndex():
    return util.getZCenterIndex()

# Initialize peaks with the best guess for height, background and sigma.
def initializePeaks(peaks, image, background, sigma, zvalue):
    c_peaks = numpy.ascontiguousarray(peaks)
    c_image = numpy.ascontiguousarray(image)
    c_background = numpy.ascontiguousarray(background)
    util.initializePeaks(c_peaks,
                         c_image,
                         c_background,
                         sigma,
                         zvalue,
                         c_peaks.shape[0],
                         c_image.shape[1])
    return c_peaks
    
# Merge new peaks with current peak list
def mergeNewPeaks(cur_peaks, new_peaks, radius, neighborhood):
    n_peak_par = getNResultsPar()
    c_cur_peaks = numpy.ascontiguousarray(cur_peaks)
    c_new_peaks = numpy.ascontiguousarray(new_peaks)
    num_cur_peaks = c_cur_peaks.size/n_peak_par
    num_new_peaks = c_new_peaks.size/n_peak_par
    c_out_peaks = numpy.ascontiguousarray(numpy.zeros(n_peak_par*(num_cur_peaks+num_new_peaks)))
    num_added = util.mergeNewPeaks(c_cur_peaks,
                                   c_new_peaks,
                                   c_out_peaks,
                                   radius,
                                   neighborhood,
                                   num_cur_peaks,
                                   num_new_peaks)
    total_peaks = n_peak_par*(num_cur_peaks+num_added)
    c_out_peaks.resize(total_peaks)
    return numpy.reshape(c_out_peaks, (-1, n_peak_par))

# Calculate the distance to the nearest peaks in (x2, y2) from (x1, y1).
def peakToPeakDist(x1, y1, x2, y2):
    c_x1 = numpy.ascontiguousarray(x1).astype(numpy.float64)
    c_y1 = numpy.ascontiguousarray(y1).astype(numpy.float64)
    c_x2 = numpy.ascontiguousarray(x2).astype(numpy.float64)
    c_y2 = numpy.ascontiguousarray(y2).astype(numpy.float64)
    n_x1 = x1.size
    n_x2 = x2.size
    c_dist = numpy.ascontiguousarray(numpy.zeros(n_x1))
    util.peakToPeakDist(c_x1, c_y1, c_x2, c_y2, c_dist, n_x1, n_x2)
    return c_dist

# Calculate the distance to the nearest peaks in (x2, y2) from (x1, y1).
def peakToPeakIndex(x1, y1, x2, y2):
    c_x1 = numpy.ascontiguousarray(x1).astype(numpy.float64)
    c_y1 = numpy.ascontiguousarray(y1).astype(numpy.float64)
    c_x2 = numpy.ascontiguousarray(x2).astype(numpy.float64)
    c_y2 = numpy.ascontiguousarray(y2).astype(numpy.float64)
    n_x1 = x1.size
    n_x2 = x2.size
    c_index = numpy.ascontiguousarray(numpy.zeros(n_x1)).astype(numpy.int32)
    util.peakToPeakIndex(c_x1, c_y1, c_x2, c_y2, c_index, n_x1, n_x2)
    return c_index

# Remove peaks that are too close to a bright neighbor from the list
def removeClosePeaks(peaks, radius, neighborhood):
    n_peak_par = getNResultsPar()
    c_in_peaks = numpy.ascontiguousarray(peaks)
    c_out_peaks = numpy.ascontiguousarray(numpy.zeros(n_peak_par*(c_in_peaks.size)))
    num_c_in_peaks = c_in_peaks.size/n_peak_par
    num_left = util.removeClosePeaks(c_in_peaks,
                                     c_out_peaks,
                                     radius,
                                     neighborhood,
                                     num_c_in_peaks)
    total_peaks = n_peak_par*num_left
    c_out_peaks.resize(total_peaks)
    return numpy.reshape(c_out_peaks, (-1, n_peak_par))

# Remove peaks that are too close to their neighbors
def removeNeighbors(peaks, radius):
    n_peak_par = getNResultsPar()
    c_in_peaks = numpy.ascontiguousarray(peaks)
    c_out_peaks = numpy.ascontiguousarray(numpy.zeros(n_peak_par*(c_in_peaks.size)))
    num_c_in_peaks = c_in_peaks.size/n_peak_par
    num_left = util.removeNeighbors(c_in_peaks,
                                    c_out_peaks,
                                    radius,
                                    num_c_in_peaks)
    total_peaks = n_peak_par*num_left
    c_out_peaks.resize(total_peaks)
    return numpy.reshape(c_out_peaks, (-1, n_peak_par))

# Smooth an image by convolving with a gaussian
def smoothImage(image):
    temp = numpy.ascontiguousarray(numpy.copy(image))
    util.smoothImage(temp,
                     temp.shape[0])
    return temp


if __name__ == "__main__":
    if 0:
        to_test = [["Background index:", getBackgroundIndex],
                   ["Error index:", getErrorIndex],
                   ["Height index:", getHeightIndex],
                   ["Number peak params:", getNPeakPar],
                   ["Number results params:", getNResultsPar],
                   ["Status index:", getStatusIndex],
                   ["X center index:", getXCenterIndex],
                   ["X width index:", getXWidthIndex],
                   ["Y center index:", getYCenterIndex],
                   ["Y width index:", getYWidthIndex],
                   ["Z center index:", getZCenterIndex]]
        for test in to_test:
            print test[0], test[1]()
    
    if 0:
        peaks = numpy.array([[11.0, 10.0, 1.0, 10.0, 1.0, 0.0, 0.0, 1, 0.0],
                             [10.0, 11.0, 1.0, 10.0, 1.0, 0.0, 0.0, 1, 0.0],
                             [11.0, 12.0, 1.0, 10.0, 1.0, 0.0, 0.0, 1, 0.0]])
        print peaks.shape
        print removeClosePeaks(peaks, 1.5, 1.5)

    if 0:
        peaks = numpy.array([[11.0, 10.0, 1.0, 10.0, 1.0, 0.0, 0.0, 1, 0.0],
                             [11.0, 12.0, 1.0, 10.0, 1.0, 0.0, 0.0, 1, 0.0],
                             [11.0, 14.0, 1.0, 10.0, 1.0, 0.0, 0.0, 1, 0.0]])
        new_peaks = numpy.array([[11.0, 12.0, 1.0, 12.0, 1.0, 0.0, 0.0, 0, 0.0]])        
        print mergeNewPeaks(peaks, new_peaks, 2.5, 4.0)

    if 0:
        size = 50
        image = numpy.zeros((size,size))
        taken = numpy.zeros((size,size), dtype=numpy.int32)
        image[20,20] = 1.0
        image[24,20] = 1.0
        image[20,24] = 1.0
        #image[size/4,size/2] = 1.0
        #image[size/2,size/4] = 1.0
        #image[size/2,size/2] = 1.0
        for i in range(3):
            print numpy.max(taken)
            [peaks, taken] = findLocalMaxima(image, 
                                             taken, 
                                             0.5,
                                             5.0,
                                             0.1,
                                             1.0, 
                                             10)
            print peaks

    if 0:
        x1 = numpy.arange(4)
        y1 = numpy.arange(4)
        x2 = numpy.arange(5)+0.1
        y2 = numpy.arange(5)
        x1[2:] += 1.0
        y1[2:] += 1.0
        print x1
        print y1
        print x2
        print y2
        print peakToPeakDist(x1, y1, x2, y2)

    if 1:
        x1 = numpy.array([1,2,3,4,5])
        y1 = numpy.array([1,2,3,4,5])
        x2 = numpy.array([3,4,1,2])
        y2 = numpy.array([3,4,1,2])
        print peakToPeakIndex(x1, y1, x2, y2)


#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
