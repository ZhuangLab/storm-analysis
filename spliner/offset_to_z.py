#!/usr/bin/python
#
# Convert a storm-control acquired movie .off file into a format
# appropriate for use by measure_psf_beads.py.
#
# Hazen 1/16
#

import numpy
import sys

if (len(sys.argv) != 3):
    print "usage: offset_to_z <.off file, input> <z file, output>"
    exit()

data = numpy.loadtxt(sys.argv[1], skiprows = 1)

offset = data[:,1]
stagez = data[:,3]

# First pass fit.
f1 = numpy.polyfit(offset, stagez, 1)
p1 = numpy.poly1d(f1)
x1 = numpy.array([numpy.min(offset) - 1, numpy.max(offset) + 1])
y1 = numpy.array(p1(x1))

# Identify outliers.
diff = stagez - p1(offset)
diff_std = numpy.std(diff)
mask = (numpy.abs(diff) < (3.0 * diff_std))

# Second pass fit.
f2 = numpy.polyfit(offset[mask], stagez[mask], 1)
p2 = numpy.poly1d(f2)
y2 = numpy.array(p2(x1))

# Testing
if 1:
    import matplotlib
    import matplotlib.pyplot as pyplot

    pyplot.figure()
    #pyplot.scatter(offset[mask], diff[mask], color = 'blue')
    #pyplot.scatter(offset[~mask], diff[~mask], color = 'red')
    pyplot.scatter(offset, stagez)
    pyplot.plot(x1, y1, color = 'red')
    pyplot.plot(x1, y2, color = 'black')
    pyplot.show()
    
