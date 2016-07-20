#!/usr/bin/python
#
# Given a spline file create an object that generate PSFs
# at requested z values.
#
# Hazen 01/16
#

import pickle
import numpy

import spline3D

class SplineToPSF(object):

    def __init__(self, spline_file):
        spline_data = pickle.load(open(spline_file))
        self.zmin = spline_data["zmin"]
        self.zmax = spline_data["zmax"]
        self.spline = spline3D.Spline3D(spline_data["spline"], spline_data["coeff"])
        self.spline_size = self.spline.getSize()

    def getPSF(self, z_value, shape = None, up_sample = 1, normalize = True):

        # Calculate PSF at requested z value.
        scaled_z = self.getScaledZ(z_value)

        psf_size = up_sample * (self.spline_size - 1)/2
        psf = numpy.zeros((psf_size, psf_size))
        for x in range(psf_size):
            for y in range(psf_size):
                psf[y,x] = self.spline.f(scaled_z,
                                         float(2*y)/float(up_sample) + 1.0,
                                         float(2*x)/float(up_sample) + 1.0)

        # Draw into a larger image if requested.
        if shape is not None:
            im_size_x = shape[0] * up_sample
            im_size_y = shape[1] * up_sample

            start_x = im_size_x/2 - psf_size/2
            start_y = im_size_y/2 - psf_size/2
            end_x = start_x + psf_size
            end_y = start_y + psf_size

            temp = numpy.zeros((im_size_x, im_size_y))
            temp[start_x:end_x,start_y:end_y] = psf

            psf = temp

        # Normalize if requested.
        if normalize:
            psf = psf/numpy.sum(psf)
            
        return psf

    def getScaledZ(self, z_value):
        return float(self.spline_size) * (z_value - self.zmin) / (self.zmax - self.zmin)

    def getSize(self):
        return self.spline_size

    def getZMin(self):
        return self.zmin

    def getZMax(self):
        return self.zmax
    

if (__name__ == "__main__"):
    import sys
    import sa_library.daxwriter as daxwriter

    if (len(sys.argv) != 3):
        print "usage: <spline (input)> <dax (output)>"
        exit()

    stp = SplineToPSF(sys.argv[1])
    size = (stp.getSize() - 1)/2
    dax_data = daxwriter.DaxWriter(sys.argv[2], size, size)
    for z in [-500.0, -250.0, 0.0, 250.0, 500.0]:
        psf = stp.getPSF(z)
        dax_data.addFrame(1000.0 * psf + 100.0)

    dax_data.close()
    
    
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
