#!/usr/bin/env python
"""
Creates a localization file with a random distribution of emitters.

Hazen 09/17
"""
import argparse
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py

parser = argparse.ArgumentParser(description = "Create a uniform random distribution of emitters for simulations.")

parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                    help = "The name of the HDF5 file to save the emitter locations, etc.")
parser.add_argument('--density', dest='density', type=float, required=True,
                    help = "Localizations per pixel squared.")
parser.add_argument('--margin', dest='margin', type=int, required=False, default = 10,
                    help = "The margin in pixels around the edge of the image, default is 10.")
parser.add_argument('--sx', dest='sx', type=int, required=False, default = 256,
                    help = "The image size in Y in pixels, default is 256.")
parser.add_argument('--sy', dest='sy', type=int, required=False, default = 256,
                    help = "The image size in Y in pixels, default is 256.")
parser.add_argument('--zrange', dest='zrange', type=float, required=False, default = 0.0,
                    help = "Range for z values in microns, -zrange to zrange, default is 0.0.")

args = parser.parse_args()

# Set random number generator seed.
numpy.random.seed(0)

# Calculate the size of area covered by the localizations.
size_x = args.sx - 2*args.margin
size_y = args.sy - 2*args.margin

# Calculate number of localizations.
n_locs = int(round(size_x*size_y*args.density))

# Create localizations.
peaks = {}
peaks["x"] = args.margin + size_y * numpy.random.uniform(size = n_locs)
peaks["y"] = args.margin + size_x * numpy.random.uniform(size = n_locs)
peaks["z"] = -args.zrange + 2.0 * args.zrange * numpy.random.uniform(size = n_locs)
peaks["xsigma"] = 1.5*numpy.ones(peaks["x"].size)
peaks["ysigma"] = 1.5*numpy.ones(peaks["y"].size)

# Save localizations.
with saH5Py.SAH5Py(args.hdf5, is_existing = False, overwrite = True) as h5:
    h5.setMovieProperties(1,1,1,"")
    h5.addLocalizations(peaks, 0)
