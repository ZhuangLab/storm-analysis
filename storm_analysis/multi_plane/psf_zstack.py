#!/usr/bin/env python
"""
Given a movie and list of locations (the output of
multi_plane.psf_localizations), generate an average
z stack.

Hazen 05/17
"""

import numpy
import scipy
import scipy.ndimage
import tifffile

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.readinsight3 as readinsight3


def psfZStack(movie_name, i3_filename, zstack_name, aoi_size = 8):

    # Load movie
    movie_data = datareader.inferReader(movie_name)
    movie_len = movie_data.filmSize()[2]
    
    # Load localizations
    i3_data = readinsight3.loadI3File(i3_filename)
    x = i3_data["x"]
    y = i3_data["y"]

    z_stack = numpy.zeros((4*aoi_size, 4*aoi_size, movie_len))

    for i in range(movie_len):
        if ((i%50) == 0):
            print("Processing frame", i)
            
        frame = movie_data.loadAFrame(i)

        for j in range(x.size):
            xf = x[j]
            yf = y[j]
            xi = int(xf)
            yi = int(yf)

            im_slice = frame[xi - aoi_size:xi + aoi_size,
                             yi - aoi_size:yi + aoi_size]

            im_slice_up = scipy.ndimage.interpolation.zoom(im_slice, 2.0)
            im_slice_up = scipy.ndimage.interpolation.shift(im_slice_up, (-2.0*(xf-xi), -2.0*(yf-yi)), mode='nearest')

            z_stack[:,:,i] += im_slice_up

    z_stack = z_stack.astype(numpy.float32)
    
    # Save z_stack
    with tifffile.TiffWriter(zstack_name) as tf:
        for i in range(movie_len):
            tf.save(z_stack[:,:,i])

            
if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Average AOIs together into a z-stack for PSF measurement.')

    parser.add_argument('--movie', dest='movie', type=str, required=False,
                        help = "The name of the movie, can be .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations psf file.")
    parser.add_argument('--zstack', dest='zstack', type=str, required=True,
                        help = "The name of the file to save the zstack in (.tif).")
    parser.add_argument('--aoi_size', dest='aoi_size', type=int, required=False, default=8,
                        help = "The size of the area of interest around the bead in pixels. The default is 8.")

    args = parser.parse_args()
    
    psfZStack(args.movie,
              args.mlist,
              args.zstack,
              aoi_size = args.aoi_size)
