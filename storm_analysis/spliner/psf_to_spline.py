#!/usr/bin/python
#
# Takes a 3D numpy array (as created by measure_psf.py) and
# outputs a 3D numpy array that can be used as a spline for
# 3D spline fitting.
#
# Size is the size of the spline in pixels. This should be
# a multiple of 2.
#
# Hazen 01/14
#

import pickle
import numpy
import sys

import spline1D
import spline2D
import spline3D

if (len(sys.argv)!=4):
    print("usage: psf_to_spline <psf file, input> <spline file, output> <size, input>")
    exit()

psf_data = pickle.load(open(sys.argv[1]))
np_psf = psf_data["psf"]
s_size = int(sys.argv[3])
spline = False
start = np_psf.shape[1]/2.0 - s_size - 0.5


# 2D spline
if (len(np_psf.shape) == 2):
    print("Generating 2D spline.")
    s_size = 2*s_size

    np_spline = numpy.zeros((s_size, s_size))
    #np_psf = np_psf/numpy.max(np_psf)
    xy_spline = spline2D.Spline2D(np_psf)

    x = start
    for i in range(s_size):
        y = start
        for j in range(s_size):
            np_spline[j,i] = xy_spline.f(y,x)
            
            y += 1.0
        x += 1.0

    print("Calculating spline coefficients.")
    spline = spline2D.Spline2D(np_spline)

    if 1:
        import sa_library.daxwriter as daxwriter
        daxwriter.singleFrameDax("spline.dax", 1000.0*np_spline + 100)


# 3D spline
else:
    print("Generating 3D spline.")
    s_size = 2*s_size

    np_spline = numpy.zeros((s_size, s_size, s_size))
    xy_splines = []

    print("Generating XY splines.")
    for i in range(np_psf.shape[0]):
        xy_splines.append(spline2D.Spline2D(np_psf[i,:,:]))

    print("Generating fitting spline.")
    x = start
    for i in range(s_size):
        y = start
        for j in range(s_size):

            zvals = numpy.zeros(np_psf.shape[0])
            for k in range(np_psf.shape[0]):
                zvals[k] = xy_splines[k].f(y,x)
            z_spline = spline1D.Spline1D(zvals)

            max_z = float(np_psf.shape[0]) - 1.0
            inc = max_z/(float(s_size)-1.0)
            for k in range(s_size):
                z = float(k)*inc
                if (z > max_z):
                    z = max_z
                np_spline[k,j,i] = z_spline.f(z)

            y += 1.0
        x += 1.0

    print("Calculating spline coefficients.")
    spline = spline3D.Spline3D(np_spline)

    if 1:
        import sa_library.daxwriter as daxwriter
        dxw = daxwriter.DaxWriter("spline.dax", np_spline.shape[1], np_spline.shape[2])
        for i in range(s_size):
            dxw.addFrame(1000.0*np_spline[i,:,:] + 100)
        dxw.close()

del psf_data["psf"]
psf_data["spline"] = np_spline
psf_data["coeff"] = spline.getCoeff()
pickle.dump(psf_data, open(sys.argv[2], "w"))

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
