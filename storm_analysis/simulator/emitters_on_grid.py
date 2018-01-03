#!/usr/bin/env python
"""
Creates lists of molecules on a grid with a +-0.5 pixel 
random offset.

Hazen 12/16
"""
import argparse
import numpy
import random

import storm_analysis.sa_library.sa_h5py as saH5Py

parser = argparse.ArgumentParser(description = "Create a grid of emitters for simulations.")

parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                    help = "The name of the HDF5 file to save the emitter locations, etc.")
parser.add_argument('--nx', dest='nx', type=int, required=True,
                    help = "The grid size in X.")
parser.add_argument('--ny', dest='ny', type=int, required=True,
                    help = "The grid size in Y.")
parser.add_argument('--spacing', dest='spacing', type=float, required=True,
                    help = "The grid spacing in pixels.")
parser.add_argument('--zrange', dest='zrange', type=float, required=False, default = 0.0,
                    help = "Range for z values in nm, -zrange to zrange")
parser.add_argument('--zoffset', dest='zoffset', type=float, required=False, default = 0.0,
                    help = "Offset for z values in nm")

args = parser.parse_args()

random.seed(0)

nx = args.nx
ny = args.ny
spacing = args.spacing
z_range = args.zrange

if (nx*ny > 1):
    curz = -z_range
    z_inc = 2.0 * z_range/(nx*ny - 1)
else:
    curz = 0.0
    z_inc = 0.0

peaks = {"id" : numpy.zeros(nx*ny, dtype = numpy.int32),
         "x" : numpy.zeros(nx*ny),
         "y" : numpy.zeros(nx*ny),
         "z" : numpy.zeros(nx*ny),
         "xsigma" : numpy.ones(nx*ny),
         "ysigma" : numpy.ones(nx*ny)}

curx = spacing
for i in range(nx):
    cury = spacing
    for j in range(ny):
        k = i*ny+j
        peaks['x'][k] = cury + random.random() - 0.5
        peaks['y'][k] = curx + random.random() - 0.5
        peaks['z'][k] = curz + args.zoffset

        # Record emitter id in the 'id' field.
        peaks['id'][k] = k
        
        cury += spacing
        curz += z_inc
    curx += spacing

with saH5Py.SAH5Py(args.hdf5, is_existing = False, overwrite = True) as h5:
    h5.setMovieProperties(1,1,1,"")
    h5.addLocalizations(peaks, 0)
