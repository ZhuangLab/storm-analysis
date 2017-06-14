#!/usr/bin/env python
"""
Giving a mapping file (from multi_plane.mapper), and a 
molecule list, generate molecule lists to use for the
PSF extraction step.

Hazen 05/17
"""

import numpy
import os
import pickle

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.ia_utilities_c as util_c
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3


def psfLocalizations(i3_filename, mapping_filename, frame = 1, aoi_size = 8, movie_filename = None):

    # Load localizations.
    i3_reader = readinsight3.I3Reader(i3_filename)

    # Load mapping.
    mappings = {}
    if os.path.exists(mapping_filename):
        with open(mapping_filename, 'rb') as fp:
            mappings = pickle.load(fp)
    else:
        print("Mapping file not found, single channel data?")

    # Try and determine movie frame size.
    i3_metadata = readinsight3.loadI3Metadata(i3_filename)
    if i3_metadata is None:
        if movie_filename is None:
            raise Exception("I3 metadata not found and movie filename is not specified.")
        else:
            movie_fp = datareader.inferReader(movie_filename)
            [movie_y, movie_x] = movie_fp.filmSize()[:2]
    else:
        movie_data = i3_metadata.find("movie")

        # FIXME: These may be transposed?
        movie_x = int(movie_data.find("movie_x").text)
        movie_y = int(movie_data.find("movie_y").text)
    
    # Load localizations in the requested frame.
    locs = i3_reader.getMoleculesInFrame(frame)
    print("Loaded", locs.size, "localizations.")

    # Remove localizations that are too close to each other.
    in_locs = numpy.zeros((locs["x"].size, util_c.getNPeakPar()))
    in_locs[:,util_c.getXCenterIndex()] = locs["x"]
    in_locs[:,util_c.getYCenterIndex()] = locs["y"]

    out_locs = util_c.removeNeighbors(in_locs, 2 * aoi_size)

    xf = out_locs[:,util_c.getXCenterIndex()]
    yf = out_locs[:,util_c.getYCenterIndex()]

    #
    # Remove localizations that are too close to the edge or
    # outside of the image in any of the channels.
    #
    is_good = numpy.ones(xf.size, dtype = numpy.bool)
    for i in range(xf.size):

        # Check in Channel 0.
        if (xf[i] < aoi_size) or (xf[i] + aoi_size >= movie_x):
            is_good[i] = False
            continue
        
        if (yf[i] < aoi_size) or (yf[i] + aoi_size >= movie_y):
            is_good[i] = False
            continue

        # Check other channels.
        for key in mappings:
            if not is_good[i]:
                break
            
            coeffs = mappings[key]
            [ch1, ch2, axis] = key.split("_")
            if (ch1 == "0"):

                if (axis == "x"):
                    xm = coeffs[0] + coeffs[1]*xf[i] + coeffs[2]*yf[i]
                    if (xm < aoi_size) or (xm + aoi_size >= movie_x):
                        is_good[i] = False
                        break

                elif (axis == "y"):
                    ym = coeffs[0] + coeffs[1]*xf[i] + coeffs[2]*yf[i]
                    if (ym < aoi_size) or (ym + aoi_size >= movie_y):
                        is_good[i] = False
                        break

    #
    # Save localizations for each channel.
    #
    gx = xf[is_good]
    gy = yf[is_good]

    basename = os.path.splitext(i3_filename)[0]
    with writeinsight3.I3Writer(basename + "_c0_psf.bin") as w3:
        w3.addMoleculesWithXY(gx, gy)
    
    index = 1
    while ("0_" + str(index) + "_x" in mappings):
        cx = mappings["0_" + str(index) + "_x"]
        cy = mappings["0_" + str(index) + "_y"]
        #cx = mappings[str(index) + "_0" + "_x"]
        #cy = mappings[str(index) + "_0" + "_y"]
        xm = cx[0] + cx[1] * gx + cx[2] * gy
        ym = cy[0] + cy[1] * gx + cy[2] * gy

        with writeinsight3.I3Writer(basename + "_c" + str(index) + "_psf.bin") as w3:
            w3.addMoleculesWithXY(xm, ym)

        index += 1

    #
    # Print localizations that were kept.
    #
    print(gx.size, "localizations were kept:")
    for i in range(gx.size):
        print("ch0: {0:.2f} {1:.2f}".format(gx[i], gy[i]))
        index = 1
        while ("0_" + str(index) + "_x" in mappings):
            cx = mappings["0_" + str(index) + "_x"]
            cy = mappings["0_" + str(index) + "_y"]
            xm = cx[0] + cx[1] * gx[i] + cx[2] * gy[i]
            ym = cy[0] + cy[1] * gx[i] + cy[2] * gy[i]
            print("ch" + str(index) + ": {0:.2f} {1:.2f}".format(xm, ym))
            index += 1
        print("")
    print("")


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Determine localizations to use for PSF measurement.')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations file. This is a binary file in Insight3 format.")
    parser.add_argument('--map', dest='mapping', type=str, required=True,
                        help = "The name of the mapping file. This is the output of multi_plane.mapper.")
    parser.add_argument('--frame', dest='frame', type=int, required=False, default=1,
                        help = "The frame in .bin file to get the localizations from. The default is 1.")
    parser.add_argument('--aoi_size', dest='aoi_size', type=int, required=False, default=8,
                        help = "The size of the area of interest around the bead in pixels. The default is 8.")
    parser.add_argument('--movie', dest='movie', type=str, required=False,
                        help = "The name of the movie, can be .dax, .tiff or .spe format.")

    args = parser.parse_args()
    
    psfLocalizations(args.mlist,
                     args.mapping,
                     frame = args.frame,
                     aoi_size = args.aoi_size,
                     movie_filename = args.movie)
    
