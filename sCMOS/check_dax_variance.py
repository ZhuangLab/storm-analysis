#!/usr/bin/python
#
# Calculate the pixel variance of a dax movie. This is useful
# for verifying that you are using the right camera
# calibration, assuming that your dax movie does not have a
# lot of real signal overwhelming the read out noise of
# the camera.
#
# Hazen 10/13
#

import numpy
import sys

import sa_library.datareader as datareader

if (len(sys.argv) != 3):
    print "usage: <input_dax> <variance>"
    exit()

cam_offset = 100
max_frames = 1000

# Open the input file.
in_file = datareader.inferReader(sys.argv[1])
[w, h, l] = in_file.filmSize()

if (l > max_frames):
    l = max_frames

# Calculate x and xx.
mean = numpy.zeros((w,h), dtype = numpy.int64)
var = numpy.zeros((w,h), dtype = numpy.int64)

for i in range(l):
    if ((i%10)==0):
        print "Processing frame", i

    aframe = in_file.loadAFrame(i)

    aframe = aframe.astype(numpy.int64)
    aframe -= cam_offset

    mean += aframe
    var += aframe * aframe

# Calculate mean and variance
mean = mean.astype(numpy.float64)
var = var.astype(numpy.float64)
mean = mean/float(l)
var = var/float(l) - mean*mean

numpy.save(sys.argv[2], [mean, var])

#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
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
