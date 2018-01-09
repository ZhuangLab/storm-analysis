#!/usr/bin/env python
"""
Creates lists of molecules in clusters

Hazen 09/17
"""
import argparse
import numpy
import random

import storm_analysis.sa_library.sa_h5py as saH5Py

parser = argparse.ArgumentParser(description = "Create emitters in (possibly overlapping) clusters.")

parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                    help = "The name of the HDF5 file to save the emitter locations, etc.")
parser.add_argument('--ncl', dest='ncl', type=int, required=True,
                    help = "The number of clusters.")
parser.add_argument('--nlocs', dest='nlocs', type=int, required=True,
                    help = "The number of localizations per cluster.")
parser.add_argument('--dev', dest='dev', type=float, required=True,
                    help = "Cluster standard deviation in pixels.")
parser.add_argument('--sx', dest='sx', type=int, required=False, default=256,
                    help = "Image x size in pixels, default is 256.")
parser.add_argument('--sy', dest='sy', type=int, required=False, default=256,
                    help = "Image y size in pixels, default is 256.")
parser.add_argument('--z_start', dest='z_start', type=int, required=False, default=-500,
                    help = "Starting value for z position, default is -500nm.")
parser.add_argument('--z_stop', dest='z_stop', type=int, required=False, default=500,
                    help = "Stopping value for z position, default is 500nm.")

args = parser.parse_args()

# First, create a list of cluster centers.
cl_centers = []
while (len(cl_centers) < args.ncl):
    cx = random.uniform(0.0, args.sx)
    cy = random.uniform(0.0, args.sy)
    cz = random.uniform(args.z_start, args.z_stop)

    # Don't keep the cluster if it is too close to the edge of the image.
    if (cx < 2.0) or (cx > (args.sx - 2.0)):
        continue
    if (cy < 2.0) or (cy > (args.sy - 2.0)):
        continue

    cl_centers.append([cx, cy, cz])

# Next, create localizations for each cluster.
xp = None
yp = None
zp = None
for clc in cl_centers:

    if xp is None:
        xp = numpy.random.normal(scale = args.dev, size = args.nlocs) + clc[0]
        yp = numpy.random.normal(scale = args.dev, size = args.nlocs) + clc[1]

        # Z is in nm, we'll assume a 100nm pixel size.
        zp = numpy.random.normal(scale = args.dev * 100.0, size = args.nlocs) + clc[2]
    else:
        xp = numpy.append(xp, numpy.random.normal(scale = args.dev, size = args.nlocs) + clc[0])
        yp = numpy.append(yp, numpy.random.normal(scale = args.dev, size = args.nlocs) + clc[1])
        zp = numpy.append(zp, numpy.random.normal(scale = args.dev * 100.0, size = args.nlocs) + clc[2])

# Create a molecule list structure & save it.
peaks = {}
peaks["x"] = xp
peaks["y"] = yp
peaks["z"] = zp
peaks["xsigma"] = 1.5*numpy.ones(xp.size)
peaks["ysigma"] = 1.5*numpy.ones(yp.size)

# Save localizations.
saH5Py.saveLocalizations(args.hdf5, peaks)
