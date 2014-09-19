#!/usr/bin/python
#
# Measure the PSF for L1-Homotopy fitting.
#
# Hazen 10/12
#

import numpy
import scipy
import scipy.ndimage
import sys

import sa_library.ia_utilities_c as util_c
import sa_library.datareader as datareader
import sa_library.readinsight3 as readinsight3

if (len(sys.argv)!=4):
    print "usage: homotopy_psf <dax_file, input> <bin_file, input> <npy_file, output>"
    exit()

# Minimum number of peaks to calculate the PSF from.
min_peaks = 500

# Half width of the aoi size in pixels.
aoi_size = 8

# Load dax file and corresponding molecule list file.
dax_data = datareader.inferReader(sys.argv[1])
i3_data = readinsight3.loadI3File(sys.argv[2])

# Go through the frames identifying good peaks and adding them
# to the average psf
average_psf = numpy.zeros((4*aoi_size,4*aoi_size))
curf = 1
peaks_used = 0
total = 0.0
[dax_x, dax_y, dax_l] = dax_data.filmSize()
while (curf < dax_l) and (peaks_used < min_peaks):

    # Select localizations in current frame & not near the edges.
    mask = (i3_data['fr'] == curf) & (i3_data['x'] > aoi_size) & (i3_data['x'] < (dax_y - aoi_size - 1)) & (i3_data['y'] > aoi_size) & (i3_data['y'] < (dax_x - aoi_size - 1))
    xr = i3_data['x'][mask]
    yr = i3_data['y'][mask]
    ht = i3_data['h'][mask]

    # Remove localizations that are too close to each other.
    in_peaks = numpy.zeros((xr.size,util_c.getNResultsPar()))
    in_peaks[:,util_c.getXCenterIndex()] = xr
    in_peaks[:,util_c.getYCenterIndex()] = yr
    in_peaks[:,util_c.getHeightIndex()] = ht

    out_peaks = util_c.removeNeighbors(in_peaks, aoi_size)

    print curf, in_peaks.shape, out_peaks.shape

    # Use remaining localizations to calculate spline.
    image = dax_data.loadAFrame(curf-1).astype(numpy.float64)

    xr = out_peaks[:,util_c.getXCenterIndex()]
    yr = out_peaks[:,util_c.getYCenterIndex()]
    ht = out_peaks[:,util_c.getHeightIndex()]

    for i in range(xr.size):
        xf = xr[i]
        yf = yr[i]
        xi = int(xf)
        yi = int(yf)

        # get localization image
        mat = image[xi-aoi_size:xi+aoi_size,
                    yi-aoi_size:yi+aoi_size]

        # re-center image
        psf = scipy.ndimage.interpolation.shift(mat,(-(xf-xi),-(yf-yi)),mode='nearest')

        # zoom in by 2x
        psf = scipy.ndimage.interpolation.zoom(psf,2.0)

        # add to average psf accumulator
        average_psf += psf
        total += ht[i]

        peaks_used += 1
        
    curf += 1

average_psf = average_psf/total

# force psf to be zero (on average) at the boundaries.
if 1:
    edge = numpy.concatenate((average_psf[0,:],
                              average_psf[-1,:],
                              average_psf[:,0],
                              average_psf[:,-1]))
    average_psf -= numpy.mean(edge)

# save PSF (in image form).
if 1:
    import sa_library.daxwriter as daxwriter
    daxwriter.singleFrameDax("psf.dax", 1000.0*average_psf+100)

# save PSF (in numpy form).
numpy.save(sys.argv[3],average_psf)


#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
