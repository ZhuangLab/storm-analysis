#!/usr/bin/python
#
# Converts a dax file into "calibration" format, i.e. a format
# that can be read by camera_calibration.py for the purpose
# of camera calibration.
#
# Note that this is more for testing. Using "calibrate.py" in
# the STORM-Control project is more efficient as you don't
# have to save a massive movie first.
#
# Hazen 10/13
#

import numpy
import sys

import sa_library.datareader as datareader

if (len(sys.argv) != 3):
    print "usage: <input_dax> <calib>"
    exit()

cam_offset = 100
is_dax = True

# Open the input file.
in_file = datareader.inferReader(sys.argv[1])
[w, h, l] = in_file.filmSize()

# Calculate x & xx.
mean = numpy.zeros(w*h, dtype = numpy.int64)
var = numpy.zeros(w*h, dtype = numpy.int64)

for i in range(l):
    aframe = in_file.loadAFrame(i)

    # Undo standard dax transpose.
    if is_dax:
        aframe = numpy.transpose(aframe)

    aframe = aframe.astype(numpy.int64).flatten()
    aframe -= cam_offset

    mean += aframe
    var += aframe * aframe

# Save the results.
numpy.save(sys.argv[2], [numpy.array([l]), mean, var])

mean_mean = numpy.mean(mean)/float(l)
print "mean of mean:", mean_mean
print "mean of variance:", numpy.mean(var)/float(l) - mean_mean*mean_mean

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
