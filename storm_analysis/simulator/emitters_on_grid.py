#!/usr/bin/python
#
# Creates lists of molecules on a grid.
#
# Hazen 12/16
#

import random
import sys

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.writeinsight3 as writeinsight3


if (len(sys.argv) != 5):
    print("usage: <bin_file> <nx> <ny> <spacing>")
    exit()

random.seed(0)

nx = int(sys.argv[2])
ny = int(sys.argv[3])
spacing = float(sys.argv[4])

i3data = i3dtype.createDefaultI3Data(nx * ny)

curx = spacing
for i in range(nx):
    cury = spacing
    for j in range(ny):
        i3data['x'][i*ny+j] = curx + random.random() - 0.5
        i3data['y'][i*ny+j] = cury + random.random() - 0.5
        cury += spacing
    curx += spacing

with writeinsight3.I3Writer(sys.argv[1]) as i3w:
    i3w.addMolecules(i3data)
