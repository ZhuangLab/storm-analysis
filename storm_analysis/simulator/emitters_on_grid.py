#!/usr/bin/env python
"""
Creates lists of molecules on a grid with a +-0.5 pixel 
random offset.

Hazen 12/16
"""

import numpy
import random

import storm_analysis.sa_library.sa_h5py as saH5Py


def emittersOnGrid(h5_name, nx, ny, sigma, spacing, zrange, zoffset, seed = 0):

    if seed is not None:
        random.seed(seed)

    if (nx*ny > 1):
        curz = -zrange
        z_inc = 2.0 * zrange/(nx*ny - 1)
    else:
        curz = 0.0
        z_inc = 0.0

    peaks = {"id" : numpy.zeros(nx*ny, dtype = numpy.int32),
             "x" : numpy.zeros(nx*ny),
             "y" : numpy.zeros(nx*ny),
             "z" : numpy.zeros(nx*ny),
             "xsigma" : sigma * numpy.ones(nx*ny),
             "ysigma" : sigma * numpy.ones(nx*ny)}

    curx = spacing
    for i in range(nx):
        cury = spacing
        for j in range(ny):
            k = i*ny+j
            peaks['x'][k] = curx + random.random() - 0.5
            peaks['y'][k] = cury + random.random() - 0.5
            peaks['z'][k] = curz + zoffset
            
            # Record emitter id in the 'id' field.
            peaks['id'][k] = k
        
            cury += spacing
            curz += z_inc
        curx += spacing

    saH5Py.saveLocalizations(h5_name, peaks)


if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = "Create a grid of emitters for simulations.")
    
    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the HDF5 file to save the emitter locations, etc.")
    parser.add_argument('--nx', dest='nx', type=int, required=True,
                        help = "The grid size in X.")
    parser.add_argument('--ny', dest='ny', type=int, required=True,
                        help = "The grid size in Y.")
    parser.add_argument('--sigma', dest='sigma', type=float, required=False, default = 1.5,
                        help = "PSF sigma in pixels.")
    parser.add_argument('--spacing', dest='spacing', type=float, required=True,
                        help = "The grid spacing in pixels.")
    parser.add_argument('--zrange', dest='zrange', type=float, required=False, default = 0.0,
                        help = "Range for z values in microns, -zrange to zrange")
    parser.add_argument('--zoffset', dest='zoffset', type=float, required=False, default = 0.0,
                        help = "Offset for z values in microns")

    args = parser.parse_args()

    emittersOnGrid(args.hdf5, args.nx, args.ny, args.sigma, args.spacing, args.zrange, args.zoffset)
    
