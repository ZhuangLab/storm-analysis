#!/usr/bin/env python
"""
Convert a storm-control acquired movie .off file into a format
appropriate for use by measure_psf_beads.py. This makes the
following assumptions:

1. The first frame is at z = 0.
2. The drift in z is neglible, so the stagez position is accurate.

Hazen 1/16
"""
import numpy
import sys

def offsetToZ(offset_file, dz = 0.0):
    """
    offset_file - The name of the text file with the z offsets.
    dz - Additional offset (in nanometers) to apply to the z offset.
    """
    data = numpy.loadtxt(offset_file, skiprows = 1)

    stagez = data[:,3]

    z_offset = 1000.0 * (stagez - stagez[0]) + dz

    # This selects the center portion of the data where the stage is
    # actually scanning through z values.
    #
    start = numpy.argmin(stagez)
    stop = numpy.argmax(stagez)
    valid = numpy.ones(z_offset.size)
    valid[:start] = 0
    valid[stop:] = 0

    return (numpy.transpose(numpy.vstack((valid, z_offset))))
    
    
if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = 'Convert storm-control offset file to a Spliner friendly form.')

    parser.add_argument('--z_offsets', dest='z_offsets', type=str, required=True,
                        help = "The name of the file to save the z offsets in.")
    parser.add_argument('--sc_file', dest='sc_file', type=str, required=True,
                        help = "The name of storm-control format .off text file.")
    parser.add_argument('--deltaz', dest='deltaz', type=float, required=False, default = 0.0,
                        help = "Z offset (in nanometers) to add to the values in the .off file.")

    args = parser.parse_args()

    zdata = offsetToZ(args.sc_file, dz = args.deltaz)    
    numpy.savetxt(args.z_offsets, zdata)
