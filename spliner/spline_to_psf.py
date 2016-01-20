#!/usr/bin/python
#
# Given a spline file create an object that generate PSFs
# at requested z values.
#
# Hazen 01/16
#

import pickle
import numpy
import sys

import spline3D

class SplineToPSF(object):

    def __init__(self, spline_file):
        spline_data = pickle.load(spline_file)
        self.zmin = spline_data["zmin"]
        self.zmax = spline_data["zmax"]
        self.spline = spline3D.Spline3D(spline_data["spline"], spline_data["coeffs"])
        self.spline_size = self.spline.getSize()

    def getPSF(self, z_value):
        scaled_z = float(self.spline_size) * (self.zmax - self.zmin) * (z_value - self.zmin)
        psf = numpy.zeros((self.spline_size, self.spline_size))
        for x in range(self.spline_size):
            for y in range(self.spline_size):
                psf[x,y] = self.spline.f(scaled_z, y, z)
        return psf
        

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
