#!/usr/bin/python
#
# Utility functions that are used for drift correction.
#
# Hazen 10/14
#

import numpy
import scipy

# Save the x,y and z drift data to a file.
def saveDriftData(filename, fdx, fdy, fdz):
    frames = numpy.arange(fdx.size) + 1
    numpy.savetxt(filename,
                  numpy.column_stack((frames,
                                      -fdx, 
                                      -fdy, 
                                      fdz)),
                  fmt = "%d\t%.3f\t%.3f\t%.3f")

# Interpolate drift data to the length of the film.
def interpolateData(xvals, yvals, film_l):

    # Create spline for interpolation.
    sp = scipy.interpolate.interp1d(xvals, yvals, kind = "cubic")

    # interpolate.
    final_drift = numpy.zeros(film_l)
    i = int(xvals[0])
    while (i <= int(xvals[-1])):
        final_drift[i] = sp(i)
        i += 1

    # Linear extrapolation at the ends.
    i = int(xvals[0])
    diff = final_drift[i] - final_drift[i+1]
    cur = final_drift[i]
    while (i >= 0):
        final_drift[i] = cur
        cur += diff
        i -= 1

    i = int(xvals[-1])
    diff = final_drift[i] - final_drift[i-1]
    cur = final_drift[i]
    while (i < film_l):
        final_drift[i] = cur
        cur += diff
        i += 1

    return final_drift

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
