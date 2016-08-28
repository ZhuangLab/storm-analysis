#!/usr/bin/python
#
# Perform compressed sensing analysis on a dax file using the
# homotopy approach. Return the results in hres image format and
# as a list of object locations.
#
# Hazen 09/12
#

import numpy
import os
import sys

import homotopy_imagea_c
import sa_library.datareader as datareader
import sa_library.parameters as parameters
import sa_library.readinsight3 as readinsight3
import sa_library.writeinsight3 as writeinsight3
import setup_A_matrix

#
# Setup
#
src_directory = os.path.dirname(__file__)

if(len(sys.argv)!=5):
    print "usage: cs_analysis <dax_file> <params_file> <hres_file> <bin_file>"
    exit()

movie_data = datareader.inferReader(sys.argv[1])

#
# FIXME:
#
# This should also start at the same frame as hres in the event of a restart.
#
i3_file = writeinsight3.I3Writer(sys.argv[4])

params = parameters.Parameters(sys.argv[2])

#
# Load the a matrix and setup the homotopy image analysis class.
#
a_mat_file = params.a_matrix

print "Using A matrix file:", a_mat_file
a_mat = setup_A_matrix.loadAMatrix(a_mat_file)

image = movie_data.loadAFrame(0)
htia = homotopy_imagea_c.HomotopyIA(a_mat,
                                    params.epsilon,
                                    image.shape)

#
# This opens the file. If it already exists, then it sets the file pointer
# to the end of the file & returns the number of the last frame analyzed.
#
curf = htia.openHRDataFile(sys.argv[3])

#
# Figure out which frame to start & stop at.
#
[dax_x,dax_y,dax_l] = movie_data.filmSize()

if hasattr(params, "start_frame"):
    if (params.start_frame>=curf) and (params.start_frame<dax_l):
        curf = params.start_frame

if hasattr(params, "max_frame"):
    if (params.max_frame>0) and (params.max_frame<dax_l):
        dax_l = params.max_frame

print "Starting analysis at frame", curf

#
# Analyze the dax data.
#
total_peaks = 0
try:
    while(curf<dax_l):

        # Load image, subtract baseline & remove negative values.
        image = movie_data.loadAFrame(curf).astype(numpy.float)
        image -= params.baseline
        mask = (image < 0)
        image[mask] = 0

        # Analyze image.
        hres_image = htia.analyzeImage(image)
        peaks = htia.saveHRFrame(hres_image, curf + 1)
        [cs_x,cs_y,cs_a,cs_i] = htia.getPeaks(hres_image)
        i3_file.addMoleculesWithXYAItersFrame(cs_x, cs_y, cs_a, cs_i, curf+1)

        peaks = cs_x.size
        total_peaks += peaks
        print "Frame:", curf, peaks, total_peaks

        curf += 1

except KeyboardInterrupt:
    print "Analysis stopped."

# cleanup
htia.closeHRDataFile()
i3_file.close()

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
