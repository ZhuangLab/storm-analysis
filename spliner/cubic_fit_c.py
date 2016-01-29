#!/usr/bin/python
#
# Simple Python interface to cubic_spline.c.
#
# Hazen 01/14
#

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import sa_library.ia_utilities_c as util_c

import spline2D
import spline3D


# Load the library.
directory = os.path.dirname(__file__)
if not (directory == ""):
    directory += "/"

if (sys.platform == "win32"):
    cubic_fit = ctypes.cdll.LoadLibrary(directory + "cubic_fit.dll")
else:
    cubic_fit = ctypes.cdll.LoadLibrary(directory + "cubic_fit.so")


# C interface definition.
cubic_fit.fSpline2D.argtypes = [ctypes.c_double,
                                ctypes.c_double]

cubic_fit.fSpline2D.restype = ctypes.c_double

cubic_fit.fSpline3D.argtypes = [ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double]

cubic_fit.fSpline3D.restype = ctypes.c_double


cubic_fit.getResidual.argtypes = [ndpointer(dtype=numpy.float64)]

cubic_fit.getResults.argtypes = [ndpointer(dtype=numpy.float64)]

cubic_fit.getUnconverged.restype = ctypes.c_int

cubic_fit.getZOff.restype = ctypes.c_int

cubic_fit.getZSize.restype = ctypes.c_int

cubic_fit.initSpline2D.argtypes = [ndpointer(dtype=numpy.float64),
                                   ctypes.c_int,
                                   ctypes.c_int]

cubic_fit.initSpline3D.argtypes = [ndpointer(dtype=numpy.float64),
                                   ctypes.c_int,
                                   ctypes.c_int,
                                   ctypes.c_int]

cubic_fit.initializeMultiFit.argtypes = [ndpointer(dtype=numpy.float64),
                                         ctypes.c_double,
                                         ctypes.c_int,
                                         ctypes.c_int]

cubic_fit.newImage.argtypes = [ndpointer(dtype=numpy.float64)]

cubic_fit.newPeaks2D.argtypes = [ndpointer(dtype=numpy.float64),
                                 ctypes.c_int]

cubic_fit.newPeaks3D.argtypes = [ndpointer(dtype=numpy.float64),
                                 ctypes.c_int]


# Globals
default_tol = 1.0e-6
height_index = util_c.getHeightIndex()
n_results_par = util_c.getNResultsPar()
status_index = util_c.getStatusIndex()
z_index = util_c.getZCenterIndex()


#
# Functions
#
def fSpline2D(x, y):
    return cubic_fit.fSpline2D(x, y)

def fSpline3D(x, y, z):
    return cubic_fit.fSpline3D(x, y, z)


#
# Classes.
#
class CSplineFit():

    def __init__(self, spline_vals, scmos_data, tolerance = default_tol):
        self.initialized = False
        self.iterations = 0
        self.peaks_size = 0
        self.scmos_data = False
        self.tolerance = tolerance

        # Initialize fitter.
        if (type(scmos_data) == type(numpy.array([]))):
            self.initializeCubicFit(scmos_data)

    def cleanup(self):
        cubic_fit.splineCleanup()
        cubic_fit.multiFitCleanup()
        self.scmos_data = False

    def doFit(self, peaks, max_iterations = 200):
        self.newPeaks(peaks)
        self.iterateSpline()
        while ((self.getUnconverged() > 0) and (self.iterations < max_iterations)):
            self.iterateSpline()
       
    def freePeaks(self):
        self.peaks_size = 0
        cubic_fit.freePeaks()

    def getCoeff(self):
        return self.py_spline.getCoeff()

    def getGoodPeaks(self, min_height = 0.0, verbose = False):
        peaks = self.getPeaks()
        if (peaks.size > 0):
            mask = (peaks[:,status_index] != 2.0) & (peaks[:,height_index] > min_height)
            if verbose:
                print " ", numpy.sum(mask), "were good out of", peaks.shape[0]
            return peaks[mask,:]
        else:
            return peaks

    def getIterations(self):
        return self.iterations

    def getPeaks(self):
        peaks = numpy.zeros((self.peaks_size), dtype = numpy.float64)
        cubic_fit.getResults(peaks)
        peaks = numpy.reshape(peaks, (-1, n_results_par)) 
        return peaks

    def getResidual(self):
        residual = numpy.zeros((self.scmos_data.shape), dtype = numpy.float64)
        cubic_fit.getResidual(numpy.ascontiguousarray(residual))
        return residual

    def getSize(self):
        return self.py_spline.getSize()
        
    def getUnconverged(self):
        return cubic_fit.getUnconverged()

    def initializeCubicFit(self, scmos_data):
        self.scmos_data = scmos_data
        cubic_fit.initializeMultiFit(numpy.ascontiguousarray(scmos_data, dtype = numpy.float64),
                                     self.tolerance,
                                     scmos_data.shape[1],
                                     scmos_data.shape[0])
        self.initialized = True

    def iterateSpline(self):
        self.iterations += 1
        cubic_fit.iterateSpline()
        
    def newPeaks(self, peaks):
        self.iterations = 0
        self.peaks_size = peaks.size
        
    def newImage(self, image):
        if (not self.initialized):
            self.initializeCubicFit(numpy.zeros(image.shape))

        if (image.shape == self.scmos_data.shape):
            cubic_fit.newImage(numpy.ascontiguousarray(image, dtype = numpy.float64))
        else:
            print "image size must match scmos data size", image.shape, self.scmos_data.shape

    def rescaleZ(self, peaks, zmin, zmax):
        return peaks


class CSpline2DFit(CSplineFit):

    def __init__(self, spline_vals, coeff_vals, scmos_data, tolerance = default_tol):
        CSplineFit.__init__(self, spline_vals, scmos_data, tolerance)

        # Initialize spline.
        self.py_spline = spline2D.Spline2D(spline_vals, coeff = coeff_vals)
        cubic_fit.initSpline2D(numpy.ascontiguousarray(self.py_spline.coeff, dtype = numpy.float64),
                               self.py_spline.max_i,
                               self.py_spline.max_i)

    def newPeaks(self, peaks):
        CSplineFit.newPeaks(self, peaks)
        n_peaks = peaks.size/n_results_par
        cubic_fit.newPeaks2D(numpy.ascontiguousarray(peaks, dtype = numpy.float64),
                             n_peaks)


class CSpline3DFit(CSplineFit):

    def __init__(self, spline_vals, coeff_vals, scmos_data, tolerance = default_tol):
        CSplineFit.__init__(self, spline_vals, scmos_data, tolerance)

        # Initialize spline.
        self.py_spline = spline3D.Spline3D(spline_vals, coeff = coeff_vals)
        cubic_fit.initSpline3D(numpy.ascontiguousarray(self.py_spline.coeff, dtype = numpy.float64),
                               self.py_spline.max_i,
                               self.py_spline.max_i,
                               self.py_spline.max_i)

    def newPeaks(self, peaks):
        CSplineFit.newPeaks(self, peaks)
        n_peaks = peaks.size/n_results_par
        cubic_fit.newPeaks3D(numpy.ascontiguousarray(peaks, dtype = numpy.float64),
                             n_peaks)

    def rescaleZ(self, peaks, zmin, zmax):
        cubic_zoff = cubic_fit.getZOff()
        cubic_zrange = cubic_fit.getZSize()
        #print "rZ:", cubic_zoff, cubic_zrange
        spline_range = zmax - zmin
        inv_zscale = 1.0/float(cubic_zrange)
        peaks[:,z_index] = ((peaks[:,z_index] - cubic_zoff)*inv_zscale*spline_range) + zmin
        return peaks


# Tests.
if __name__ == "__main__":
    pass


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
