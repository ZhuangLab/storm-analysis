#!/usr/bin/python
#
# Automated XYZ drift correction for STORM movies.
#
# Hazen 1/10
#
# Modified to deal better with super huge insight3 files.
#
# Hazen 11/11
#

import numpy
import os
import scipy.signal
import sys

import sa_library.arraytoimage as arraytoimag
import sa_library.driftutilities as driftutilities
import sa_library.i3togrid as i3togrid
import sa_library.imagecorrelation as imagecorrelation


# Setup
if (len(sys.argv) < 5):
    print "usage: <bin> <drift.txt> <step> <scale> <optional - z_correct>"
    exit()

step = int(sys.argv[3])
scale = int(sys.argv[4])
i3_data = i3togrid.I3GDataLL(sys.argv[1], scale = scale)
film_l = i3_data.getFilmLength()

correct_z = True
if (len(sys.argv) > 5):
    correct_z = False

# Sub-routines.
def saveDriftData(fdx, fdy, fdz):
    driftutilities.saveDriftData(sys.argv[2], fdx, fdy, fdz)

def interpolateData(xvals, yvals):
    return driftutilities.interpolateData(xvals, yvals, film_l)

# Don't analyze films that are too short.
if ((4*step) >= film_l):
    saveDriftData(numpy.zeros(film_l),
                  numpy.zeros(film_l),
                  numpy.zeros(film_l))
    exit()

#
# Drift correction (XY and Z are all done at the same time)
#
# Note that drift corrected localizations are added back into 
# the reference image in the hopes of improving the correction
# for subsequent localizations. 
#

start = 0
i3_data.loadDataInFrames(fmin = start, fmax = start+step-1)
xymaster = i3_data.i3To2DGridAllChannelsMerged(uncorrected = True)

if correct_z:
    z_bins = 20
    xyzmaster = i3_data.i3To3DGridAllChannelsMerged(z_bins,
                                                    uncorrected = True)

index = 1
last = 0
step_step = 0
if(start>0):
    j = 0
else:
    j = step
t = [step/2]
x = [0]
y = [0]
z = [0]
old_dx = 0.0
old_dy = 0.0
old_dz = 0.0
while(j < film_l):

    # Load correct frame range.
    last = j
    if ((j + 2*step) >= film_l):
        i3_data.loadDataInFrames(fmin = j)
        step_step = 2*step
    else:
        i3_data.loadDataInFrames(fmin = j, fmax = j + step - 1)
        step_step = step

    xycurr = i3_data.i3To2DGridAllChannelsMerged(uncorrected = True)

    # Correlate to master image.
    [corr, dx, dy, xy_success] = imagecorrelation.xyOffset(xymaster,
                                                           xycurr,
                                                           i3_data.getScale(),
                                                           center = [x[index-1] * scale,
                                                                     y[index-1] * scale])

    # Update values
    if xy_success:
        old_dx = dx
        old_dy = dy
    else:
        dx = old_dx
        dy = old_dy

    dx = dx/float(scale)
    dy = dy/float(scale)

    t.append(step/2 + index * step)
    x.append(dx)
    y.append(dy)

    i3_data.applyXYDriftCorrection(dx,dy)
    if xy_success:
        # Add current to master
        xymaster += i3_data.i3To2DGridAllChannelsMerged()

    # Z correlation
    dz = old_dz
    if correct_z and xy_success:

        xyzcurr = i3_data.i3To3DGridAllChannelsMerged(z_bins,
                                                      uncorrected = True)

        # Do z correlation
        [corr, fit, dz, z_success] = imagecorrelation.zOffset(xyzmaster, xyzcurr)

        # Update Values
        if z_success:
            old_dz = dz
        else:
            dz = old_dz
            
        dz = dz * 1000.0/float(z_bins)

        if z_success:
            i3_data.applyZDriftCorrection(-dz)
            xyzmaster += i3_data.i3To3DGridAllChannelsMerged(z_bins)
    
    z.append(dz)

    print index, dx, dy, dz

    index += 1
    j += step_step

i3_data.close()

# Create numpy versions of the drift arrays.
nt = numpy.array(t)
final_driftx = interpolateData(nt, numpy.array(x))
final_drifty = interpolateData(nt, numpy.array(y))
final_driftz = interpolateData(nt, numpy.array(z))

saveDriftData(final_driftx,
              final_drifty,
              final_driftz)

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
