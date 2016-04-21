#!/usr/bin/python
#
# Calculate the error (for simulations) between the known
# peak locations and the peak locations returned by the analysis.
#
# Note: This measures the error and *not* the recall. If only 1
# in 100 peaks is found but their locations agree exactly with
# the known locations then this will be considered good.
#
# One assumption here is that neither of these files has a
# particularly large number of localizations.
#
# Hazen 01/16
#

import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import sys

import sa_library.gaussfit as gaussfit
import sa_library.readinsight3 as readinsight3
import sa_library.ia_utilities_c as utilC

if (len(sys.argv) != 3):
    print "usage: <true locations> <measured locations>"
    exit()

# For converting XY units to nanometers.
pixel_size = 160.0

truth_i3 = readinsight3.I3Reader(sys.argv[1])
measured_i3 = readinsight3.I3Reader(sys.argv[2])

all_dx = None
all_dy = None
all_dz = None
for i in range(truth_i3.getNumberFrames()):
    t_locs = truth_i3.getMoleculesInFrame(i+1)
    m_locs = measured_i3.getMoleculesInFrame(i+1, good_only = False)

    p_index = utilC.peakToPeakIndex(m_locs['xc'], m_locs['yc'], t_locs['xc'], t_locs['yc'])
    dx = numpy.zeros(m_locs.size)
    dy = numpy.zeros(m_locs.size)
    dz = numpy.zeros(m_locs.size)
    for i in range(m_locs.size):
        dx[i] = pixel_size * (m_locs['xc'][i] - t_locs['xc'][p_index[i]])
        dy[i] = pixel_size * (m_locs['yc'][i] - t_locs['yc'][p_index[i]])
        dz[i] = m_locs['zc'][i] - t_locs['zc'][p_index[i]]

    if all_dx is None:
        all_dx = dx
        all_dy = dy
        all_dz = dz
    else:
        all_dx = numpy.concatenate((all_dx, dx))
        all_dy = numpy.concatenate((all_dy, dy))
        all_dz = numpy.concatenate((all_dz, dz))

print "means and standard deviations (in nm):"
print "mean, std (dx)", numpy.mean(all_dx), numpy.std(all_dx)
print "mean, std (dy)", numpy.mean(all_dy), numpy.std(all_dy)
print "mean, std (dz)", numpy.mean(all_dz), numpy.std(all_dz)
print ""

[hist_dx, bins] = numpy.histogram(all_dx, bins = 30, range = (-100.0, 100.0))
[hist_dy, bins] = numpy.histogram(all_dy, bins = 30, range = (-100.0, 100.0))
[hist_dz, bins] = numpy.histogram(all_dz, bins = 30, range = (-100.0, 100.0))

hist_dx = hist_dx.astype(numpy.float)/numpy.sum(hist_dx)
hist_dy = hist_dy.astype(numpy.float)/numpy.sum(hist_dy)
hist_dz = hist_dz.astype(numpy.float)/numpy.sum(hist_dz)

centers = bins[:-1] + 0.5 * (bins[1] - bins[0])

print "gaussian fitting"
bin_size = bins[1] - bins[0]
[fitx, goodx] =  gaussfit.fitSymmetricGaussian1D(hist_dx)
[fity, goody] =  gaussfit.fitSymmetricGaussian1D(hist_dy)
[fitz, goodz] =  gaussfit.fitSymmetricGaussian1D(hist_dz)

print ""
print "gaussian fit to error histogram width (in nm):"
if goodx:
    print "x width", fitx[3]*bin_size
if goody:
    print "y width", fity[3]*bin_size
if goodz:
    print "z width", fitz[3]*bin_size    
    

fig = pyplot.figure()
pyplot.plot(centers, hist_dx, color = "red")
pyplot.plot(centers, hist_dy, color = "green")
pyplot.plot(centers, hist_dz, color = "blue")
pyplot.xlabel("Error in nm")
pyplot.ylabel("Density (AU)")
pyplot.show()

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
