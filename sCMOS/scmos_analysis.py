#!/usr/bin/python
#
# Perform scmos analysis on a dax file given parameters.
#
# Hazen 02/14
#

import sys

import find_peaks
import sa_library.parameters as params
import sa_utilities.std_analysis as std_analysis

# setup
if(len(sys.argv)==3):
    parameters = params.Parameters(sys.argv[2])
    mlist_file = sys.argv[1][:-4] + "_mlist.bin"
elif(len(sys.argv)==4):
    parameters = params.Parameters(sys.argv[3])
    mlist_file = sys.argv[2]
else:
    print "usage: <movie> <bin> <parameters.xml>"
    exit()

finder = find_peaks.initFindAndFit(parameters)    
std_analysis.standardAnalysis(finder,
                              sys.argv[1],
                              mlist_file,
                              parameters)

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
