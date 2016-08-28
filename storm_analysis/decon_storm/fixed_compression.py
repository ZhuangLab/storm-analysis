#!/usr/bin/python
#
# Performs iterative deconvolution with a fixed compression parameter
# and a background that is assumed to be uniform.
#
# Hazen 05/12
#

import numpy
import sys

import arraytoimage as arraytoimage
import datareader as datareader
import daxwriter as daxwriter
import mlem_c as mlem

# defaults
scale = 8
iters = 500

# user defined
input_movie = datareader.SPEReader(sys.argv[1])

[x_size, y_size, frames] = input_movie.filmSize()

if (x_size != y_size):
    print "Movies must be square.."
    exit()

output_movie = daxwriter.DaxWriter(sys.argv[2], x_size*scale, y_size*scale)
camera_offset = float(sys.argv[3])
sigma = float(sys.argv[4])
compression = float(sys.argv[5])

mlemd = mlem.Fitter(numpy.zeros((x_size,y_size)),
                    sigma,
                    scale,
                    0.0)

# process the film
for i in range(frames):
    print "Processing:", i
    
    # load image
    image = input_movie.loadAFrame(i) - camera_offset

    # deconvolve
    mlemd.newForeground(image)
    mlemd.doFitCompressedFixedBg(compression, verbose = True, max_iter = iters)
    high_res = mlemd.getForeground()
    high_res = 1000.0 * high_res/numpy.max(high_res)
    output_movie.addFrame(high_res)

output_movie.close()
