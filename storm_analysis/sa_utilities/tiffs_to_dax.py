#!/usr/bin/python
#
# Convert a mess of tiff files into a single dax file.
#
# Hazen 09/14
#

import glob
import numpy
import sys

import sa_library.daxwriter as daxwriter
import sa_library.datareader as datareader

if (len(sys.argv) != 3):
    print "usage: <dax> <tiff dir>"
    exit()

dax_file = daxwriter.DaxWriter(sys.argv[1], 0, 0)
tiff_files = sorted(glob.glob(sys.argv[2] + "*.tif"))

if (len(tiff_files) == 0):
    print "No tiff files found in '" + sys.argv[2] + "'"
    exit()

for tiff_image in tiff_files:
    print tiff_image
    data = datareader.TifReader(tiff_image).loadAFrame(0)
    if 0:
        data = data - numpy.median(data) + 2000
    dax_file.addFrame(data)

dax_file.close()

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
