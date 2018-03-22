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
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.sa_h5py as saH5Py


def psfLocalizations(h5_filename, mapping_filename, frame = 0, aoi_size = 8, min_height = 0.0):

    # Load localizations & movie size.
    with saH5Py.SAH5Py(h5_filename) as h5:
        locs = h5.getLocalizationsInFrame(frame)
        assert bool(locs), "No localizations found in frame " + str(frame)
        [movie_x, movie_y] = h5.getMovieInformation()[:2]

    # Load mapping.
    mappings = {}
    if os.path.exists(mapping_filename):
        with open(mapping_filename, 'rb') as fp:
            mappings = pickle.load(fp)
    else:
        print("Mapping file not found, single channel data?")

    # Remove localizations that are too dim.
    mask = (locs["height"] > min_height)

    locs_mask = {}
    for elt in ["x", "y"]:
        locs_mask[elt] = locs[elt][mask]
    
    # Remove localizations that are too close to each other.
    [xf, yf] = iaUtilsC.removeNeighbors(locs_mask["x"], locs_mask["y"], 2.0 * aoi_size)

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

    basename = os.path.splitext(h5_filename)[0]
    saH5Py.saveLocalizations(basename + "_c1_psf.hdf5", {"x" : gx, "y" : gy})
    
    index = 1
    while ("0_" + str(index) + "_x" in mappings):
        cx = mappings["0_" + str(index) + "_x"]
        cy = mappings["0_" + str(index) + "_y"]
        xm = cx[0] + cx[1] * gx + cx[2] * gy
        ym = cy[0] + cy[1] * gx + cy[2] * gy

        saH5Py.saveLocalizations(basename + "_c" + str(index+1) + "_psf.hdf5", {"x" : xm, "y" : ym})
        
        index += 1

    #
    # Print localizations that were kept.
    #
    print(gx.size, "localizations were kept out of", xf.size)
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
                        help = "The name of the localizations file.")
    parser.add_argument('--map', dest='mapping', type=str, required=True,
                        help = "The name of the mapping file. This is the output of multi_plane.mapper.")
    parser.add_argument('--frame', dest='frame', type=int, required=False, default=0,
                        help = "The frame in .bin file to get the localizations from. The default is 0.")
    parser.add_argument('--aoi_size', dest='aoi_size', type=int, required=False, default=8,
                        help = "The size of the area of interest around the bead in pixels. The default is 8.")
    parser.add_argument('--min_height', dest='min_height', type=float, required=False, default = 0.0,
                        help = "Minimum localization height.")

    args = parser.parse_args()
    
    psfLocalizations(args.mlist,
                     args.mapping,
                     frame = args.frame,
                     aoi_size = args.aoi_size,
                     min_height = args.min_height)
    
