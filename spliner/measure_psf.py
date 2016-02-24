#!/usr/bin/python
#
# Measure the PSF given the raw data and localization analysis. In
# theory this will be more accurate as you are using fit locations
# instead of locations that were entered by hand.
#
# 2D: The expected input is a movie file that contains images
#     of single molecules & the corresponding analysis file
#     as created by Spliner, 3D-DAOSTORM, sCMOS or Insight3.
#
# 3D: The expected input is a movie file that contains images of
#     single molecules at different z positions. This movie file
#     needs to have been analyzed with Spliner, 3D-DAOSTORM, sCMOS or 
#     Insight3.
#
# Similar to measure_psf_beads, depending on your setup you may need to change:
#  1. The z range (z_range).
#  2. The pixel size (pixel_size).
#  3. The AOI size (aoi_size). This is less important as you
#     will get to specify the final size to use when you use
#     psf_to_spline.py to create the spline to use for fitting.
#
# Hazen 01/16
#

import pickle
import numpy
import scipy
import scipy.ndimage
import sys

import sa_library.ia_utilities_c as util_c
import sa_library.datareader as datareader
import sa_library.readinsight3 as readinsight3

if (len(sys.argv)!= 5):
    print "usage: measure_psf <dax_file, input> <bin_file, input> <psf_file, output> <3d (0,1), input>"
    exit()

# Half width of the aoi size in pixels.
aoi_size = 10

# Load dax file and corresponding molecule list file.
dax_data = datareader.DaxReader(sys.argv[1])
i3_data = readinsight3.loadI3File(sys.argv[2])

# Determine whether this is 2D or 3D.
analysis_type = "3D"
if (sys.argv[4] == "0"):
    print "Measuring 2D PSF"
    analysis_type = "2D"
else:
    print "Measuring 3D PSF"

#
# Go through the frames identifying good peaks and adding them
# to the average psf. For 3D molecule z positions are rounded to 
# the nearest 50nm.
#
z_range = 500.0  # This is really more like the half range, the
                 # full range will cover +- z_range.

z_step = 50.0
z_mid = int(z_range/z_step)
max_z = 2 * z_mid + 1

average_psf = numpy.zeros((max_z,4*aoi_size,4*aoi_size))
peaks_used = 0
totals = numpy.zeros(max_z)
[dax_x, dax_y, dax_l] = dax_data.filmSize()
for curf in range(dax_l):

    # Select localizations in current frame & not near the edges.
    mask = (i3_data['fr'] == curf+1) & (i3_data['x'] > aoi_size) & (i3_data['x'] < (dax_x - aoi_size - 1)) & (i3_data['y'] > aoi_size) & (i3_data['y'] < (dax_y - aoi_size - 1))
    xr = i3_data['x'][mask]
    yr = i3_data['y'][mask]
    zr = i3_data['z'][mask]
    ht = i3_data['h'][mask]

    # Remove localizations that are too close to each other.
    in_peaks = numpy.zeros((xr.size,util_c.getNResultsPar()))
    in_peaks[:,util_c.getXCenterIndex()] = xr
    in_peaks[:,util_c.getYCenterIndex()] = yr
    in_peaks[:,util_c.getZCenterIndex()] = zr
    in_peaks[:,util_c.getHeightIndex()] = ht

    out_peaks = util_c.removeNeighbors(in_peaks, 2*aoi_size)
    #out_peaks = util_c.removeNeighbors(in_peaks, aoi_size)

    print curf, "peaks in", in_peaks.shape[0], ", peaks out", out_peaks.shape[0]

    # Use remaining localizations to calculate spline.
    image = dax_data.loadAFrame(curf).astype(numpy.float64)

    xr = out_peaks[:,util_c.getXCenterIndex()]
    yr = out_peaks[:,util_c.getYCenterIndex()]
    zr = out_peaks[:,util_c.getZCenterIndex()]
    ht = out_peaks[:,util_c.getHeightIndex()]

    for i in range(xr.size):
        xf = xr[i]
        yf = yr[i]
        zf = zr[i]
        xi = int(xf)
        yi = int(yf)
        if (analysis_type == "2D"):
            zi = 0
        else:
            zi = round(zf/z_step) + z_mid

        # check the z is in range
        if (zi > -1) and (zi < max_z):

            # get localization image
            mat = image[xi-aoi_size:xi+aoi_size,
                        yi-aoi_size:yi+aoi_size]

            # zoom in by 2x
            psf = scipy.ndimage.interpolation.zoom(mat, 2.0)

            # re-center image
            psf = scipy.ndimage.interpolation.shift(psf, (-2.0*(xf-xi), -2.0*(yf-yi)), mode='nearest')

            # add to average psf accumulator
            average_psf[zi,:,:] += psf
            totals[zi] += ht[i]

# Force PSF to be zero (on average) at the boundaries.
for i in range(max_z):
    edge = numpy.concatenate((average_psf[i,0,:],
                              average_psf[i,-1,:],
                              average_psf[i,:,0],
                              average_psf[i,:,-1]))
    average_psf[i,:,:] -= numpy.mean(edge)

# Normalize the PSF.
if (analysis_type == "2D"):
    max_z = 1

for i in range(max_z):
    print i, totals[i]
    if (totals[i] > 0.0):
        # The factor of 4.0 corrects for the 2x upsampling.
        average_psf[i,:,:] = 4.0 * average_psf[i,:,:]/numpy.sum(numpy.abs(average_psf[i,:,:]))

# Save PSF (in image form).
if 1:
    import sa_library.daxwriter as daxwriter
    dxw = daxwriter.DaxWriter("psf.dax", average_psf.shape[1], average_psf.shape[2])
    for i in range(max_z):
        dxw.addFrame(10000.0 * average_psf[i,:,:] + 100)
    dxw.close()

# Save PSF.
if (analysis_type == "2D"):
    dict = {"psf" : average_psf[0,:,:],
            "type" : "2D"}

else:
    cur_z = -z_range
    z_vals = []
    for i in range(max_z):
        z_vals.append(cur_z)
        cur_z += z_step

    dict = {"psf" : average_psf,
            "pixel_size" : 0.080, # 1/2 the camera pixel size in nm.
            "type" : "3D",
            "zmin" : -z_range,
            "zmax" : z_range,
            "zvals" : z_vals}

pickle.dump(dict, open(sys.argv[3], "w"))

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
