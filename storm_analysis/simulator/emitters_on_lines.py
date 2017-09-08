#!/usr/bin/env python
"""
Creates lists of emitters on randomly oriented lines.

Hazen 08/17
"""
import argparse
import math
import random

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.writeinsight3 as writeinsight3


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

parser.add_argument('--bin', dest='i3bin', type=str, required=True,
                    help = "The name of Insight3 format file to save the emitter locations, etc.")
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
parser.add_argument('--z_start', dest='z_start', type=int, required=False, default=-500,
                    help = "Starting value for z position, default is -500nm.")
parser.add_argument('--z_stop', dest='z_stop', type=int, required=False, default=500,
                    help = "Stopping value for z position, default is 500nm.")
        
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
i3data = i3dtype.createDefaultI3Data(args.nemitters)
i = 0
printed = False
while (i < args.nemitters):

    if((i%1000) == 0):
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
        
        i3data[i]['x'] = pnt[0]
        i3data[i]['xc'] = pnt[0]

        i3data[i]['y'] = pnt[1]
        i3data[i]['yc'] = pnt[1]

        i3data[i]['z'] = pnt[2]
        i3data[i]['zc'] = pnt[2]

        i += 1
        
# Save emitter locations.
with writeinsight3.I3Writer(args.i3bin) as i3w:
    i3w.addMolecules(i3data)

    
