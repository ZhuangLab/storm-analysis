#!/usr/bin/env python
"""
Use matplotlib to create an image of a frame with localizations overlaid.

Hazen 03/18
"""
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
import numpy

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.sa_h5py as saH5Py


def overlayImage(movie_name, locs_name, frame_number, sx = 8, sy = 8):
    """
    Create an image of a frame with the localizations overlaid.

    movie_name - The name of the movie file.
    locs_name - The name of the localizations HDF5 file.
    frame_number - Which frame to examine.
    sx - Figure x size in inches.
    sy - Figure y size in inches.
    """
    frame = datareader.inferReader(movie_name).loadAFrame(frame_number).astype(numpy.float64)
    with saH5Py.SAH5Py(locs_name) as h5:
        locs = h5.getLocalizationsInFrame(frame_number)

    frame = frame - numpy.min(frame)
    frame = frame/numpy.max(frame)
    
    fig = pyplot.figure(figsize = (sx, sy))
    ax = fig.add_subplot(1,1,1)
    ax.imshow(frame, interpolation = 'nearest', cmap = "gray")
    for i in range(locs["x"].size):
        width = 10
        height = 10
        if "xsigma" in locs:
            width = height = 5.0*locs["xsigma"][i]
        if "ysigma" in locs:
            height = 5.0*locs["ysigma"][i]
        ellipse = patches.Ellipse((locs["x"][i], locs["y"][i]), width, height, facecolor='none', edgecolor='g', linewidth = 2)
        ax.add_artist(ellipse)
        
    #ax.scatter(locs["x"], locs["y"], s = 200,
    ax.set_title("Overlay Image")

    pyplot.show()

    
def overlayImageBeads(movie_name, beads_locs_name, frame_number, sx = 8, sy = 8):
    """
    Create an image of a frame with the bead locations overlaid.

    movie_name - The name of the movie file.
    beads_locs_name - The name of the text file with the bead locations.
    frame_number - Which frame to examine.
    sx - Figure x size in inches.
    sy - Figure y size in inches.
    """
    
    frame = datareader.inferReader(movie_name).loadAFrame(frame_number).astype(numpy.float64)
    frame = frame - numpy.min(frame)
    frame = frame/numpy.max(frame)

    bead_locs = numpy.loadtxt(beads_locs_name)
    locs = {"x" : bead_locs[:,0],
            "y" : bead_locs[:,1]}
    
    fig = pyplot.figure(figsize = (sx, sy))
    ax = fig.add_subplot(1,1,1)
    ax.imshow(frame, interpolation = 'nearest', cmap = "gray")
    for i in range(locs["x"].size):
        width = 10
        height = 10
        ellipse = patches.Ellipse((locs["x"][i], locs["y"][i]), width, height, facecolor='none', edgecolor='g', linewidth = 2)
        ax.add_artist(ellipse)
        
    #ax.scatter(locs["x"], locs["y"], s = 200,
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
    
