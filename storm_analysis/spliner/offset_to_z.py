#!/usr/bin/env python
"""
Convert a storm-control acquired movie .off file into a format
appropriate for use by measure_psf_beads.py. This makes the
following assumptions:

1. The first frame is at z = 0.
2. There drift in z is neglible, so the stagez position is accurate.

Hazen 1/16
"""

import numpy
import sys

if (len(sys.argv) != 3):
    print("usage: offset_to_z <.off file, input> <z file, output>")
    exit()

data = numpy.loadtxt(sys.argv[1], skiprows = 1)

stagez = data[:,3]

z_offset = 1000.0 * (stagez - stagez[0])
valid = numpy.ones(z_offset.size)

numpy.savetxt(sys.argv[2], numpy.transpose(numpy.vstack((valid, z_offset))))

