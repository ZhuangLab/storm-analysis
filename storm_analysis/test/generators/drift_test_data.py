#!/usr/bin/env python
"""
Generate a test file for drift correction.

Hazen 12/17
"""
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py

clusters = 100
length = 5000
p_on = 0.1

dx = 0.25/length  # pixels.
dy = 0.5/length   # pixels.
dz = 0.2/length   # microns.

numpy.random.seed(0)

cx = numpy.random.uniform(low = 20.0, high = 236.0, size = clusters)
cy = numpy.random.uniform(low = 20.0, high = 236.0, size = clusters)
cz = numpy.random.uniform(low = -0.3, high = 0.3, size = clusters)

hdf5_name = "../data/test_drift.hdf5"
    
with open("../data/test_drift.txt", "w") as drift_fp:
    with saH5Py.SAH5Py(hdf5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(256, 256, length, "XYZZY")
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
            z = cz + drift_z + numpy.random.normal(scale = 0.05, size = clusters)

            peaks = {"x" : x[on], "y" : y[on], "z" : z[on]}
            h5.addLocalizations(peaks, i)
