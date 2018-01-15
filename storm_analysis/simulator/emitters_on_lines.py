#!/usr/bin/env python
"""
Creates lists of emitters on randomly oriented lines.

Hazen 08/17
"""
import argparse
import math
import numpy
import random

import storm_analysis.sa_library.sa_h5py as saH5Py

class Line(object):

    def __init__(self, x_pos = None, y_pos = None, z_pos = None, angle = None, length = None, **kwds):
        super().__init__(**kwds)
        
        self.length = length
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.z_pos = z_pos

        self.cos_t = math.cos(math.radians(angle))
        self.sin_t = math.sin(math.radians(angle))
        
    def getLocation(self, t):
        if (t > self.length):
            return None
        
        return [self.x_pos + t * self.cos_t,
                self.y_pos + t * self.sin_t,
                self.z_pos]


parser = argparse.ArgumentParser(description = "Create randomly oriented lines of emitters.")

parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                    help = "The name of the HDF5 file to save the emitter locations, etc.")
parser.add_argument('--nlines', dest='nlines', type=int, required=True,
                    help = "The number of lines.")
parser.add_argument('--nemitters', dest='nemitters', type=int, required=True,
                    help = "The number of emitters.")
parser.add_argument('--sx', dest='sx', type=int, required=False, default=256,
                    help = "Image x size in pixels, default is 256.")
parser.add_argument('--sy', dest='sy', type=int, required=False, default=256,
                    help = "Image y size in pixels, default is 256.")
parser.add_argument('--maxl', dest='maxl', type=int, required=False, default=200,
                    help = "Maximum line length in pixels, default is 200.")
parser.add_argument('--minl', dest='minl', type=int, required=False, default=20,
                    help = "Minimum line length in pixels, default is 20.")
parser.add_argument('--start_a', dest='start_a', type=int, required=False, default=0,
                    help = "Starting value for angular range, default is 0 degrees.")
parser.add_argument('--stop_a', dest='stop_a', type=int, required=False, default=360,
                    help = "Stopping value for angular range, default is 360 degrees.")
parser.add_argument('--z_start', dest='z_start', type=int, required=False, default=-0.5,
                    help = "Starting value for z position, default is -0.5um.")
parser.add_argument('--z_stop', dest='z_stop', type=int, required=False, default=0.5,
                    help = "Stopping value for z position, default is 0.5um.")
        
args = parser.parse_args()

# Create line segments.
lines = []
for i in range(args.nlines):
    lines.append(Line(x_pos = random.uniform(0.0, args.sx),
                      y_pos = random.uniform(0.0, args.sy),
                      z_pos = random.uniform(args.z_start, args.z_stop),
                      angle = random.uniform(args.start_a, args.stop_a),
                      length = random.uniform(args.minl, args.maxl)))

# Generate emitter locations.
peaks = {}
peaks["x"] = numpy.zeros(args.nemitters)
peaks["y"] = numpy.zeros(args.nemitters)
peaks["z"] = numpy.zeros(args.nemitters)
peaks["xsigma"] = 1.5*numpy.ones(args.nemitters)
peaks["ysigma"] = 1.5*numpy.ones(args.nemitters)

i = 0
printed = False
while (i < args.nemitters):

    if((i%10000) == 0):
        if not printed:
            print("Adding point", i)
            printed = True
    else:
        printed = False

    # Get a random line.
    line = random.choice(lines)

    # Get a random point on the line.
    pnt = line.getLocation(random.uniform(0.0, args.maxl))

    # Record location if valid.
    if pnt is not None:
        if (pnt[0] < 0.0) or (pnt[0] > args.sx):
            continue
        if (pnt[1] < 0.0) or (pnt[1] > args.sy):
            continue
        
        peaks['x'][i] = pnt[0]
        peaks['y'][i] = pnt[1]
        peaks['z'][i] = pnt[2]

        i += 1
        
# Save localizations.
with saH5Py.SAH5Py(args.hdf5, is_existing = False, overwrite = True) as h5:
    h5.setMovieInformation(args.sx, args.sy, 1, "")
    h5.addLocalizations(peaks, 0)


    
