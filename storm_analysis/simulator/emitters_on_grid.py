#!/usr/bin/env python
"""
Creates lists of molecules on a grid with a +-0.5 pixel 
random offset.

Hazen 12/16
"""

import argparse
import random

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.writeinsight3 as writeinsight3

parser = argparse.ArgumentParser(description = "Create a grid of emitters for simulations.")

parser.add_argument('--bin', dest='i3bin', type=str, required=True,
                    help = "The name of Insight3 format file to save the emitter locations, etc.")
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

i3data = i3dtype.createDefaultI3Data(nx * ny)

curx = spacing
for i in range(nx):
    cury = spacing
    for j in range(ny):
        k = i*ny+j
        i3data['x'][k] = curx + random.random() - 0.5
        i3data['y'][k] = cury + random.random() - 0.5
        i3data['z'][k] = curz + args.zoffset
        #i3data['x'][k] = curx - 0.5
        #i3data['y'][k] = cury - 0.5
        #i3data['z'][k] = -250.0

        i3data['xc'][k] = i3data['x'][k]
        i3data['yc'][k] = i3data['y'][k]
        i3data['zc'][k] = i3data['z'][k]

        # Record emitter id in the 'i' field.
        i3data['i'][k] = k
        
        cury += spacing
        curz += z_inc
    curx += spacing

with writeinsight3.I3Writer(args.i3bin) as i3w:
    i3w.addMolecules(i3data)
