#!/usr/bin/python
#
# Measure the 3D PSF a movie given the locations of the beads
# of interest in the movie and the z-offset of each frame of
# the movie. It is assumed that the drift over the time
# course of the movie is neglible.
#
# Depending on your setup you may need to change:
#  1. The z range (z_range).
#  2. The pixel size (pixel_size).
#  3. The AOI size (aoi_size). This is less important as you
#     will get to specify the final size to use when you use
#     psf_to_spline.py to create the spline to use for fitting.
#
# Hazen 1/16
#

import pickle
import numpy
import scipy
import scipy.ndimage
import sys

import sa_library.datareader as datareader

if (len(sys.argv)!= 5):
    print "usage: measure_psf_beads <movie_file, input> <z_file, input> <bead_file, input> <psf_file output>"
    exit()

# Half width of the aoi size in pixels.
aoi_size = 8

# Load movie file.
movie_data = datareader.inferReader(sys.argv[1])

#
# Load the z-offset information for the dax file.
#
#   This is a text file with one line per frame that contains the 
#   z-offset (in nm) for that frame. Each line is a space separated
#   valid, z_pos pair. If valid if 0 the frame will be ignored,
#   otherwise it will be used.
#
data = numpy.loadtxt(sys.argv[2])
valid = data[:,0]
z_off = data[:,1]

#
# Load the locations of the beads.
#
#   This is a text file the contains the locations of the beads that 
#   will be used to construct the PSF. Each line is a space separated 
#   x, y pair of bead locations (in pixels).
#
#   One way to create this file is to look at the bead movie with
#   visualizer.py and record the center positions of several beads.
#
data = numpy.loadtxt(sys.argv[3])
bead_x = data[:,0]
bead_y = data[:,1]

#
# Go through the frames and the bead images to the average psf. Z 
# positions are rounded to the nearest 50nm. You might need to 
# adjust z_range depending on your experiment.
#
#z_range = 1500.0
z_range = 550.0

z_step = 50.0
z_mid = int(z_range/z_step)
max_z = 2 * z_mid + 1
average_psf = numpy.zeros((max_z,4*aoi_size,4*aoi_size))
totals = numpy.zeros(max_z)
[dax_x, dax_y, dax_l] = movie_data.filmSize()
for curf in range(dax_l):

    if ((curf%50)==0):
        print "Processing frame:", curf

    if (abs(valid[curf]) < 1.0e-6):
    #    print "skipping", valid[curf]
        continue

    # Use bead localization to calculate spline.
    image = movie_data.loadAFrame(curf).astype(numpy.float64)

    # Get frame z and check that it is in range.
    zf = z_off[curf]
    zi = round(zf/z_step) + z_mid
    if (zi > -1) and (zi < max_z):

        for i in range(bead_x.size):

            xf = bead_x[i]
            yf = bead_y[i]
            xi = int(xf)
            yi = int(yf)

            # Get localization image.
            mat = image[xi-aoi_size:xi+aoi_size,
                        yi-aoi_size:yi+aoi_size]

            # Zoom in by 2x.
            psf = scipy.ndimage.interpolation.zoom(mat, 2.0)

            # Re-center image.
            psf = scipy.ndimage.interpolation.shift(psf, (-2.0*(xf-xi), -2.0*(yf-yi)), mode='nearest')

            # Add to average psf accumulator.
            average_psf[zi,:,:] += psf
            totals[zi] += 1

# Force PSF to be zero (on average) at the boundaries.
for i in range(max_z):
    edge = numpy.concatenate((average_psf[i,0,:],
                              average_psf[i,-1,:],
                              average_psf[i,:,0],
                              average_psf[i,:,-1]))
    average_psf[i,:,:] -= numpy.mean(edge)

# Normalize PSF.
psf_sum = numpy.sum(average_psf)
for i in range(max_z):
    if (totals[i] > 0.0):
        #average_psf[i,:,:] = average_psf[i,:,:]/(totals[i] * psf_sum)
        average_psf[i,:,:] = average_psf[i,:,:]/totals[i]
        #average_psf[i,:,:] = average_psf[i,:,:]/numpy.max(average_psf[i,:,:])

average_psf = average_psf/numpy.max(average_psf)

# Save PSF (in image form).
if 1:
    import sa_library.daxwriter as daxwriter
    dxw = daxwriter.DaxWriter("psf.dax", average_psf.shape[1], average_psf.shape[2])
    for i in range(max_z):
        dxw.addFrame(1000.0 * average_psf[i,:,:] + 100)
    dxw.close()

# Save PSF. 
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

pickle.dump(dict, open(sys.argv[4], "w"))

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
