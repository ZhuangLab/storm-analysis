#!/usr/bin/env python
#
# Given arrays of x, y, z and intensities, return an array
# of objects that emulates a double-helical PSF in a format
# compatible with drawgaussians.
#
# Hazen 01/16
#

import numpy

psf_type = "double-helix"

# Double helix PSF range.
z_min = -500.0
z_max = 500.0

sigma = 1.0
def PSF(x, y, z, h):
    num_objects = x.size
    objects = numpy.zeros((2 * num_objects, 5))
    for i in range(num_objects):
        # Should there be a warning if z is not between z_min and z_max?
        #angle = numpy.pi * (0.1 + 0.8 * ((z[i] - z_min)/(z_max - z_min)))
        angle = numpy.pi * 0.5 * ((z[i] - z_min)/(z_max - z_min))
        dx = 2.0 * sigma * numpy.cos(angle)
        dy = 2.0 * sigma * numpy.sin(angle)
        objects[2*i,:] = [x[i] + dx, y[i] + dy, h[i], sigma, sigma]
        objects[2*i+1,:] = [x[i] - dx, y[i] - dy, h[i], sigma, sigma]

    return objects

def PSFIntegral(z, h):
    return numpy.pi * 4.0 * h * sigma * sigma


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
