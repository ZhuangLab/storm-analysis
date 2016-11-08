#!/usr/local/bin/python
#
# Create a test .bin file for evaluating RCC correction.
#
# Hazen 10/14
#

import numpy
import sys

import storm_analysis.sa_library.writeinsight3 as writeinsight3

if (len(sys.argv) != 2):
    print("usage: <out.bin>")
    exit()

i3_out = writeinsight3.I3Writer(sys.argv[1])
n_locs = 100

def addMol(x,y,frame):
    x += numpy.random.normal(scale = 0.05, size = n_locs)
    y += numpy.random.normal(scale = 0.05, size = n_locs)
    i3_out.addMoleculesWithXYFrame(x,y,i)

x1 = 255.0 * numpy.random.uniform(size = n_locs)
y1 = 255.0 * numpy.random.uniform(size = n_locs)

x2 = 255.0 * numpy.random.uniform(size = n_locs)
y2 = 255.0 * numpy.random.uniform(size = n_locs)

dx = 0.0
dy = 0.0
step = 50
for i in range(6 * step):

    if ((i%100)==0):
        print("frame", i)
    j = i/step

    if 0:
        addMol(x1 + dx, y1 + dy, i)
        addMol(x2 + dx, y2 + dy, i)

    if 0:
        if (j!=2):
            addMol(x1 + dx, y1 + dy, i)
            addMol(x2 + dx, y2 + dy, i)

    if 1:
        if(j==1):
            addMol(x1 + dx, y1 + dy, i)
        elif(j==2):
            addMol(x2 + dx, y2 + dy, i)
        else:
            addMol(x1 + dx, y1 + dy, i)
            addMol(x2 + dx, y2 + dy, i)
    
    if 0:
        if((j%3)==0):
            addMol(x1 + dx, y1 + dy, i)
            addMol(x2 + dx, y2 + dy, i)
        elif((j%3)==1):
            addMol(x1 + dx, y1 + dy, i)
        else:
            addMol(x2 + dx, y2 + dy, i)

    dx += 0.0
    dy += 0.0

i3_out.close()

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
