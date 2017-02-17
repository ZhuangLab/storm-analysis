#!/usr/bin/env python
"""
Generate a test file for clustering analysis.

Hazen 02/17
"""

import numpy

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.writeinsight3 as writeinsight3

background = 5
clusters = 100
length = 20000
p_on = 0.1

numpy.random.seed(0)

cx = numpy.random.uniform(low = 20.0, high = 236.0, size = clusters)
cy = numpy.random.uniform(low = 20.0, high = 236.0, size = clusters)
cz = numpy.random.uniform(low = -300.0, high = 300.0, size = clusters)

with writeinsight3.I3Writer("../data/test_clustering_list.bin") as i3_fp:
    for i in range(length):
        
        if((i % 500) == 0):
            print("Creating frame", i)

        on = (numpy.random.uniform(size = clusters) < p_on)
        number_on = numpy.count_nonzero(on)

        x = cx + numpy.random.normal(scale = 0.5, size = clusters)
        y = cy + numpy.random.normal(scale = 0.5, size = clusters)
        z = cz + numpy.random.normal(scale = 50.0, size = clusters)

        # Add background
        x = numpy.concatenate((x[on], numpy.random.uniform(high = 256.0, size = background)))
        y = numpy.concatenate((y[on], numpy.random.uniform(high = 256.0, size = background)))
        z = numpy.concatenate((z[on], numpy.random.uniform(low = -500.0, high = 500.0, size = background)))
        
        i3d = i3dtype.createDefaultI3Data(number_on + background)
        i3dtype.posSet(i3d, "x", x)
        i3dtype.posSet(i3d, "y", y)
        i3dtype.posSet(i3d, "z", z)
        i3dtype.setI3Field(i3d, "fr", i+1)

        i3_fp.addMolecules(i3d)
