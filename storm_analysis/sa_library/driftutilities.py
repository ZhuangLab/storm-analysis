#!/usr/bin/env python
"""
Utility functions that are used for drift correction.

Hazen 02/17
"""

import numpy
import scipy


def saveDriftData(filename, fdx, fdy, fdz):
    """
    Save the x,y and z drift data to a file.
    """
    frames = numpy.arange(fdx.size) + 1
    numpy.savetxt(filename,
                  numpy.column_stack((frames,
                                      -fdx, 
                                      -fdy, 
                                      fdz)),
                  fmt = "%d\t%.3f\t%.3f\t%.3f")


def interpolateData(xvals, yvals, film_l):
    """
    Interpolate drift data to the length of the film.
    """

    final_drift = numpy.zeros(film_l)
    
    # Use polyfit for extrapolation at the end points.
    pe = numpy.poly1d(numpy.polyfit(xvals[0:2], yvals[0:2], 1))
    for i in range(int(xvals[0])):
        final_drift[i] = pe(i)

    pe = numpy.poly1d(numpy.polyfit(xvals[-2:], yvals[-2:], 1))
    for i in range(int(xvals[-1]), film_l):
        final_drift[i] = pe(i)        

    # Create linear spline for interpolation.
    sp = scipy.interpolate.interp1d(xvals, yvals, kind = "linear")

    # Interpolate.
    i = int(xvals[0])
    while (i <= int(xvals[-1])):
        final_drift[i] = sp(i)
        i += 1

    return final_drift

#
# The MIT License
#
# Copyright (c) 2017 Zhuang Lab, Harvard University
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
