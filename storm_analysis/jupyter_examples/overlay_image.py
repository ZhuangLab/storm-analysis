#!/usr/bin/env python
"""
Use matplotlib to create an image of a frame with localizations overlaid.

Hazen 03/18
"""
import matplotlib
import matplotlib.pyplot as pyplot
import numpy

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.sa_h5py as saH5Py

def overlayImage(movie_name, locs_name, frame_number, sx = 8, sy = 8):
    frame = datareader.inferReader(movie_name).loadAFrame(frame_number).astype(numpy.float64)
    with saH5Py.SAH5Py(locs_name) as h5:
        locs = h5.getLocalizationsInFrame(frame_number)

    frame = frame - numpy.min(frame)
    frame = frame/numpy.max(frame)
    
    fig = pyplot.figure(figsize = (sx, sy))
    ax = fig.add_subplot(1,1,1)
    ax.imshow(frame, interpolation = 'nearest', cmap = "gray")
    ax.scatter(locs["x"], locs["y"], s = 200, facecolors='none', edgecolors='g', linewidth = 2)
    ax.set_title("Overlay Image")

    pyplot.show()


if (__name__ == "__main__"):
    
    import argparse

    parser = argparse.ArgumentParser(description = 'Overlay localizations on an image')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the localizations file.")
    parser.add_argument('--frame', dest='frame', type=int, required=True,
                        help = "Which frame.")

    args = parser.parse_args()

    overlayImage(args.movie, args.hdf5, args.frame)
    
