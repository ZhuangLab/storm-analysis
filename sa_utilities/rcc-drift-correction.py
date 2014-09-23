#!/usr/bin/python
#
# A Python implementation of the drift algorithm described in this reference:
#
# "Localization events-based sample drift correction for localization microscopy with redundant cross-correlation algorithm", 
# Wang et al. Optics Express, 30 June 2014, Vol. 22, No. 13, DOI:10.1364/OE.22.015982.
#
# This uses the above algorithm for XY correction, then falls back to old
# approach for the Z correction.
#
# Hazen 09/14
#

import numpy
import os
import pickle
import scipy.signal
import sys

import sa_library.arraytoimage as arraytoimage
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
max_err = 0.2 * scale

if film_l < step:
    sys.exit()

correct_z = True
if (len(sys.argv) > 5):
    correct_z = False

"""
#film_l = 2000

# Compute offsets between all pairs of sub images.
endpost = film_l - step/2
old_start1 = -1
start1 = 0
end1 = start1 + step
start2 = start1
end2 = start2 + step
i = 0
j = 0
centers = [(end1 - start1)/2 + start1]
pairs = []
while (start1 < endpost):

    if (start2 > endpost):
        i += 1
        j = i
        start1 += step
        end1 = start1 + step
        start2 = start1
        end2 = start2 + step
        if (end1 > endpost):
            end1 = film_l
        if (end2 > endpost):
            end2 = film_l
        if (start1 < endpost):
            centers.append((end1 - start1)/2 + start1)

    if (start1 > endpost):
        continue

    if not (start1 == start2):
        if (old_start1 != start1):
            i3_data.loadDataInFrames(fmin = start1, fmax = end1)
            sub1 = i3_data.i3To2DGridAllChannelsMerged(uncorrected = True)
            old_start1 = start1

        i3_data.loadDataInFrames(fmin = start2, fmax = end2)
        sub2 = i3_data.i3To2DGridAllChannelsMerged(uncorrected = True)

        [corr, dx, dy, success] = imagecorrelation.xyOffset(sub1,
                                                            sub2,
                                                            i3_data.getScale())

        print "offset between frame ranges ", start1, "-" , end1 , " and ", start2, "-", end2
    #print i,j
        print " -> ", dx, dy, success
        print ""

        pairs.append([i, j, dx, dy, success])

    j += 1
    start2 += step
    end2 = start2 + step
    if (end2 > endpost):
        end2 = film_l

print "--"
print film_l

with open("test.dat", "w") as fp:
    pickle.dump([centers, pairs], fp)

"""

with open("test.dat") as fp:
    [centers, pairs] = pickle.load(fp)

# Prepare rij_x, rij_y, A matrix.
rij_x = numpy.zeros(len(pairs), dtype = numpy.float32)
rij_y = numpy.zeros(len(pairs), dtype = numpy.float32)
A = numpy.zeros((len(pairs),len(centers)-1), dtype = numpy.float32)
for i, pair in enumerate(pairs):
    rij_x[i] = pair[2]
    rij_y[i] = pair[3]
    A[i,pair[0]:pair[1]] = 1.0

# Calculate drift (pass1). 
# dx and dy contain the optimal offset between sub image i and sub image i+1 in x/y.
pinv_A = numpy.linalg.pinv(A)
dx = numpy.dot(pinv_A, rij_x)
dy = numpy.dot(pinv_A, rij_y)

# Calculate errors.
err_x = numpy.dot(A, dx) - rij_x
err_y = numpy.dot(A, dy) - rij_y

err_d = numpy.sqrt(err_x * err_x + err_y * err_y)
arg_sort_err = numpy.argsort(err_d)

# Remove bad values.

#max_err = 0.02
j = len(arg_sort_err) - 1
while (err_d[arg_sort_err[j]] > max_err):
    index = arg_sort_err[j]
    delA = numpy.delete(A, index, 0)
    if (numpy.linalg.matrix_rank(A) == (len(centers)-1)):
        print "removing", index
        A = delA
        rij_x = numpy.delete(rij_x, index, 0)
        rij_y = numpy.delete(rij_y, index, 0)
        arg_sort_err[(arg_sort_err > index)] -= 1
    j -= 1

# Calculate drift (pass2). 
pinv_A = numpy.linalg.pinv(A)
dx = numpy.dot(pinv_A, rij_x)
dy = numpy.dot(pinv_A, rij_y)

# Integrate to get final drift.
driftx = numpy.zeros((dx.size))
drifty = numpy.zeros((dy.size))
for i in range(dx.size):
    driftx[i] = numpy.sum(dx[0:i])
    drifty[i] = numpy.sum(dy[0:i])

if 1:
    for i in range(driftx.size):
        print i, driftx[i], drifty[i]


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
