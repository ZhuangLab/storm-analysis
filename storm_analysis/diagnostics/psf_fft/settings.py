#!/usr/bin/env python
"""
Settings to use in PSF FFT simulations.

Hazen 10/17
"""

camera_gain = 1.0
camera_offset = 100.0
camera_variance = 2.5
iterations = 20
margin = 1
n_frames = 10
nx = 14
ny = 9
#nx = 1
#ny = 1
#peak_locations = "peaks.hdf5"
photons = [[20, 500], [20, 1000]]
#photons = [[20, 4000]]
pixel_size = 100.0
psf_size = 30
psf_z_range = 0.6
test_z_range = 0.3
test_z_offset = 0.0
tolerance = 0.3
x_size = 300
y_size = 200
z_step = 0.2
#zmn = []
zmn = [[1.3, 2, 2]]
#zmn = [[1.3, -1, 3], [1.3, -2, 2]]
