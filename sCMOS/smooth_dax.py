#!/usr/bin/python
#
# Outputs a smoothed dax file, this can be useful for
# figuring out what thresholds to use for peak finding.
#
# Hazen 10/13
#

import numpy
import sys

import sa_library.datareader as datareader
import sa_library.daxwriter as daxwriter

import scmos_utilities_c

if (len(sys.argv) != 6):
    print "usage: <input_dax> <output_dax> <calib> <sigma> <frames>"
    exit()

# Open the input file.
in_file = datareader.inferReader(sys.argv[1])
f_len = in_file.filmSize()[2]
if (int(sys.argv[5]) > 0) and (int(sys.argv[5]) < f_len):
    f_len = int(sys.argv[5])

# Open the output file.
out_file = daxwriter.DaxWriter(sys.argv[2], 0, 0)

# Load camera calibration (sliced as appropriate 
# for the ROI) and create the smoother class.
[offset, variance, gain] = numpy.load(sys.argv[3])
smoother = scmos_utilities_c.Smoother(offset, variance, gain)

# Load images, smooth & output.
sigma_psf = int(round(float(sys.argv[4])))
for i in range(f_len):
    print "Smoothing frame", i
    in_image = in_file.loadAFrame(i)
    sm_image = smoother.smoothImage(in_image, sigma_psf) + 100.0
    out_file.addFrame(sm_image)

out_file.close()

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
