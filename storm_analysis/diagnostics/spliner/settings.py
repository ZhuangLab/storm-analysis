#!/usr/bin/env python
"""
Settings to use in Spliner simulations.

Hazen 09/17
"""

camera_gain = 1.0
camera_offset = 100.0
camera_variance = 2.5
fit_error_model = "MLE"
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
smooth_psf = False
smooth_psf_sigma = 0.5
spline_size = 10
spline_z_range = 0.75
test_z_range = 0.3
test_z_offset = 0.0
tolerance = 0.3
x_size = 300
y_size = 200
#zmn = []
zmn = [[1.3, 2, 2]]
#zmn = [[1.3, -1, 3], [1.3, -2, 2]]

