#!/usr/bin/python
#
# Generate simulated. The basic idea is that you provide a list of
# localizations in Insight3 .bin format and these are used to
# generate a series of images using the following steps:
#
# Initialization:
#   1. locs = readinsight3.loadI3File(mlist_file)
#   2. bg = background.Background(locs)
#   3. camera = camera.Camera()
#   4. pp = photophysics.Photophysics(locs, intensity)
#   5. psf = psf.PSF()
#
#
# Generation:
#   1. image = numpy.zeros((256, 256))
#   2. image += bg.getBackground()
#   3. cur_locs = pp.getEmitters()
#   4. image += psf.getPSFs(cur_locs)
#   5. image = camera.readImage(image)
#   6. saveimage(image)
#   7. saveloc(cur_locs, frame_number)
#
# Hazen 11/16
#

import numpy
import sys

import storm_analysis.sa_library.daxwriter as daxwriter
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf


if (len(sys.argv) != 5):
    print("usage: <dax_file> <bin_file> <frames> <intensity>")
    exit()

#
# Initialization.
#

# Image size.
x_size = 256
y_size = 256

dax_data = daxwriter.DaxWriter(sys.argv[1], x_size, y_size)
i3_data_in = readinsight3.loadI3File(sys.argv[2])
i3_data_out = writeinsight3.I3Writer(sys.argv[1][:-4] + "_olist.bin")
n_frames = int(sys.argv[3])
intensity = float(sys.argv[4])

# Change these as needed for your simulation.
bg = background.UniformBackground(x_size, y_size, i3_data_in)
cam = camera.Ideal(x_size, y_size, 100.0)
pp = photophysics.AlwaysOn(x_size, y_size, i3_data_in, intensity)
psf = psf.GaussianPSF(x_size, y_size, 160.0)

#
# Generate frames
#
for i in range(n_frames):
    print("Generating frame:", i)

    # Generate the new image.
    image = numpy.zeros((x_size, y_size))
    image += bg.getBackground(i)

    cur_i3 = pp.getEmitters(i)
    image += psf.getPSFs(cur_i3)

    image = cam.readImage(image)

    # Save the image.
    dax_data.addFrame(image)

    # Save the molecule locations.
    cur_i3['fr'] = i + 1
    i3_data_out.addMolecules(cur_i3)

dax_data.close()
i3_data_out.close()

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
