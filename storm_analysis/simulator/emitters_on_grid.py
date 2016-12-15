#!/usr/bin/python
#
# Creates lists of molecules on a grid.
#
# Hazen 12/16
#

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

args = parser.parse_args()

random.seed(0)

nx = args.nx
ny = args.ny
spacing = args.spacing

i3data = i3dtype.createDefaultI3Data(nx * ny)

curx = spacing
for i in range(nx):
    cury = spacing
    for j in range(ny):
        i3data['x'][i*ny+j] = curx + random.random() - 0.5
        i3data['y'][i*ny+j] = cury + random.random() - 0.5
        i3data['xc'][i*ny+j] = i3data['x'][i*ny+j]
        i3data['yc'][i*ny+j] = i3data['y'][i*ny+j]

        # Record emitter id in the 'i' field.
        i3data['i'][i*ny+j] = i*ny+j
        
        cury += spacing
    curx += spacing

with writeinsight3.I3Writer(args.i3bin) as i3w:
    i3w.addMolecules(i3data)
