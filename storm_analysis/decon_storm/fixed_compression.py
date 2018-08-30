#!/usr/bin/python
#
# Performs iterative deconvolution with a fixed compression parameter
# and a background that is assumed to be uniform.
#
# Hazen 05/12
#

import numpy
import sys

import storm_analysis.sa_library.arraytoimage as arraytoimage
import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.datawriter as datawriter

import mlem_c as mlem

# defaults
scale = 8
iters = 500

# user defined
input_movie = datareader.SPEReader(sys.argv[1])

[x_size, y_size, frames] = input_movie.filmSize()

if (x_size != y_size):
    print("Movies must be square..")
    exit()

output_movie = datawriter.DaxWriter(sys.argv[2])
camera_offset = float(sys.argv[3])
sigma = float(sys.argv[4])
compression = float(sys.argv[5])

mlemd = mlem.Fitter(numpy.zeros((x_size,y_size)),
                    sigma,
                    scale,
                    0.0)

# process the film
for i in range(frames):
    print("Processing:", i)
    
    # load image
    image = input_movie.loadAFrame(i) - camera_offset

    # deconvolve
    mlemd.newForeground(image)
    mlemd.doFitCompressedFixedBg(compression, verbose = True, max_iter = iters)
    high_res = mlemd.getForeground()
    high_res = 1000.0 * high_res/numpy.max(high_res)
    output_movie.addFrame(high_res)

output_movie.close()
