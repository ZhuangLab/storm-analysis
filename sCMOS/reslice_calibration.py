#!/usr/bin/python
#
# Reslices a camera calibration file. This is used to
# create a calibration that matches the camera ROI of
# the acquired data.
#
# Hazen 10/13
#

import numpy
import sys

if (len(sys.argv) != 7):
    print "usage: <input_calib> <output_calib> <x_start> <y_start> <x_width> <y_width>"
    exit()

# This may need to be changed to match the camera calibration data dimensions.
calib_xsize = 2048
calib_ysize = 2048
#calib_xsize = 200
#calib_ysize = 300

# Load the data & reshape.
[offset, variance, gain] = numpy.load(sys.argv[1])
offset = numpy.reshape(offset, (calib_ysize, calib_xsize))
variance = numpy.reshape(variance, (calib_ysize, calib_xsize))
gain = numpy.reshape(gain, (calib_ysize, calib_xsize))

# Slice out the ROI.
x_start = int(sys.argv[3])
y_start = int(sys.argv[4])
x_stop = x_start + int(sys.argv[5])
y_stop = y_start + int(sys.argv[6])
rs_offset = offset[y_start:y_stop,x_start:x_stop]
rs_variance = variance[y_start:y_stop,x_start:x_stop]
rs_gain = gain[y_start:y_stop,x_start:x_stop]

# Transpose (since we are also transposing the images for historical reasons).
rs_offset = numpy.transpose(rs_offset)
rs_variance = numpy.transpose(rs_variance)
rs_gain = numpy.transpose(rs_gain)    

# Save sliced calibration.
numpy.save(sys.argv[2], [rs_offset, rs_variance, rs_gain])

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
