#!/usr/bin/env python
#
# Given arrays of x, y, z and intensities, return an array
# of objects that emulates an astigmatic PSF in a format
# compatible with drawgaussians.
#
# Hazen 01/16
#

import numpy

import sa_library.multi_fit_c as multi_fit_c

psf_type = "astigmatic"

# Parameters for astigmatic PSF.
wx_params = numpy.array([2.0,  0.150, 0.40, 0.0, 0.0])
wy_params = numpy.array([2.0, -0.150, 0.40, 0.0, 0.0])

def PSF(x, y, z, h):
    num_objects = x.size
    objects = numpy.zeros((num_objects, 5))
    for i in range(num_objects):
        [sx, sy] = multi_fit_c.calcSxSy(wx_params, wy_params, z[i] * 0.001)
        objects[i,:] = [x[i], y[i], h[i], sx, sy]

    return objects

    
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
