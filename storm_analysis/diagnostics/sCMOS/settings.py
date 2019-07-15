#!/usr/bin/env python
"""
Settings to use in sCMOS simulations.

Hazen 09/17
"""

camera_gain = 1.0
camera_offset = 100.0
camera_variance = 2.5
fit_error_model = "MLE"
iterations = 20
model = "2d"
n_frames = 100
nx = 14
ny = 9
#peak_locations = "peaks.txt"
photons = [[20, 500], [20, 1000]]
pixel_size = 100.0
roi_size = 9
sensitivity_correction = 1
threshold = 6.0
tolerance = 0.3
verbosity = 1
x_size = 300
y_size = 200
