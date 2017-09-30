#!/usr/bin/env python
"""
Creates a localization file with a random distribution of emitters.

Hazen 09/17
"""

import argparse
import numpy

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.writeinsight3 as writeinsight3

parser = argparse.ArgumentParser(description = "Create a uniform random distribution of emitters for simulations.")

parser.add_argument('--bin', dest='i3bin', type=str, required=True,
                    help = "The name of Insight3 format file to save the emitter locations, etc.")
parser.add_argument('--density', dest='density', type=float, required=True,
                    help = "Localizations per pixel squared.")
parser.add_argument('--margin', dest='margin', type=int, required=False, default = 10,
                    help = "The margin in pixels around the edge of the image, default is 10.")
parser.add_argument('--sx', dest='sx', type=int, required=False, default = 256,
                    help = "The image size in Y in pixels, default is 256.")
parser.add_argument('--sy', dest='sy', type=int, required=False, default = 256,
                    help = "The image size in Y in pixels, default is 256.")
parser.add_argument('--zrange', dest='zrange', type=float, required=False, default = 0.0,
                    help = "Range for z values in nm, -zrange to zrange, default is 0.0.")

args = parser.parse_args()

# Set random number generator seed.
numpy.random.seed(0)

# Calculate the size of area covered by the localizations.
size_x = args.sx - 2*args.margin
size_y = args.sy - 2*args.margin

# Calculate number of localizations.
n_locs = int(round(size_x*size_y*args.density))

# Create localization structure.
i3data = i3dtype.createDefaultI3Data(n_locs)

# Create localizations.
i3dtype.posSet(i3data, "x", args.margin + size_x * numpy.random.uniform(size = n_locs))
i3dtype.posSet(i3data, "y", args.margin + size_y * numpy.random.uniform(size = n_locs))
i3dtype.posSet(i3data, "z", -args.zrange + 2.0 * args.zrange * numpy.random.uniform(size = n_locs))

# Save localizations.
with writeinsight3.I3Writer(args.i3bin) as i3w:
    i3w.addMolecules(i3data)
