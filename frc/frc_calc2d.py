#!/usr/bin/env python
#
# Calculate 2D FRC following Nieuwenhuizen, Nature Methods, 2013
#
# Note that this calculates the uncorrected FRC.
#
# Hazen 10/14
#

import math
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import sys

import frc_c
import sa_library.arraytoimage as arraytoimage
import sa_library.i3togrid as i3togrid
import sa_library.readinsight3 as readinsight3

pixel_size = 160.0
storm_scale = 8

if (len(sys.argv) != 3):
    print "usage: <in_list.bin> <results.txt>"
    exit()

# Load the data.
i3_grid = i3togrid.I3GData(sys.argv[1], scale = storm_scale)

# Split the data (approximately) in half & generate 2D histograms.
print "Searching for mid-point"

# For simulations the .dax file might not actually have as many 
# frames as the molecule list so use a hack to get the number of
# frames in the molecule list.
max_f = int(numpy.max(i3_grid.i3data['fr'])) + 1
locs = round(numpy.sum(i3_grid.i3To2DGridAllChannelsMerged(fmax = max_f)))

start = 0
end = max_f
half_locs = locs/2
while ((end - start) > 1):
    mid = (end - start)/2 + start
    print "  ", start, mid, end
    grid1 = i3_grid.i3To2DGridAllChannelsMerged(fmin = 0, fmax = mid)
    if (numpy.sum(grid1) < half_locs):
        start = mid
    else:
        end = mid

print " mid-point:", end    
grid1 = i3_grid.i3To2DGridAllChannelsMerged(fmin = 0, fmax = end)
grid2 = i3_grid.i3To2DGridAllChannelsMerged(fmin = end, fmax = max_f)

# Compute FFT
print "Calculating"
grid1_fft = numpy.fft.fftshift(numpy.fft.fft2(grid1))
grid2_fft = numpy.fft.fftshift(numpy.fft.fft2(grid2))

grid1_fft_sqr = grid1_fft * numpy.conj(grid1_fft)
grid2_fft_sqr = grid2_fft * numpy.conj(grid2_fft)
grid1_grid2 = grid1_fft * numpy.conj(grid2_fft)

if 1:
    arraytoimage.singleColorImage(numpy.abs(grid1_fft), "grid1")
    arraytoimage.singleColorImage(numpy.abs(grid2_fft), "grid2")

[frc, frc_counts] = frc_c.frc(grid1_fft, grid2_fft)

# Plot results
for i in range(frc.size):
    if (frc_counts[i] > 0):
        frc[i] = frc[i]/float(frc_counts[i])
    else:
        frc[i] = 0.0

xvals = numpy.arange(frc.size)
xvals = xvals/(float(grid1_fft.shape[0]) * pixel_size * (1.0/float(storm_scale)))
frc = numpy.real(frc)

fp = open(sys.argv[2], "w")
for i in range(xvals.size):
    fp.write(str(xvals[i]) + "," + str(frc[i]) + "\n")
fp.close()

fig = pyplot.figure()
ax = fig.add_subplot(111)
ax.scatter(xvals, frc)
pyplot.ylim([-0.2,1.2])
pyplot.show()

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
