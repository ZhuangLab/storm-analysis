#!/usr/bin/env python
"""
Generate a test file for drift correction.

Hazen 02/17
"""

import numpy

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.writeinsight3 as writeinsight3

clusters = 100
length = 20000
p_on = 0.1

dx = 1.0/length   # pixels.
dy = 2.0/length   # pixels.
dz = 200.0/length # nanometers.

numpy.random.seed(0)

cx = numpy.random.uniform(low = 20.0, high = 236.0, size = clusters)
cy = numpy.random.uniform(low = 20.0, high = 236.0, size = clusters)
cz = numpy.random.uniform(low = -300.0, high = 300.0, size = clusters)

with open("test_drift.txt", "w") as drift_fp:
    with writeinsight3.I3Writer("test_drift_mlist.bin") as i3_fp:
        for i in range(length):
        
            if((i % 500) == 0):
                print("Creating frame", i)

            drift_x = i * dx
            drift_y = i * dy
            drift_z = i * dz

            drift_fp.write("{0:0d}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n".format(i+1, drift_x, drift_y, drift_z))
            
            on = (numpy.random.uniform(size = clusters) < p_on)
            number_on = numpy.count_nonzero(on)

            x = cx + drift_x + numpy.random.normal(scale = 0.5, size = clusters)
            y = cy + drift_y + numpy.random.normal(scale = 0.5, size = clusters)
            z = cz + drift_z + numpy.random.normal(scale = 50.0, size = clusters)
            
            i3d = i3dtype.createDefaultI3Data(number_on)
            i3dtype.posSet(i3d, "x", x[on])
            i3dtype.posSet(i3d, "y", y[on])
            i3dtype.posSet(i3d, "z", z[on])

            i3_fp.addMolecules(i3d)
