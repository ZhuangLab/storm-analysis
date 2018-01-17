#!/usr/bin/env python
"""
Settings to use in Multicolor simulations.

Note: Background photons are per plane, total photons are divided across all the planes.

Hazen 01/18
"""
import numpy

camera_gain = 1.0
camera_offset = 100.0
camera_variance = 1.0

dx = 2.0
dy = 1.0

independent_heights = 1
iterations = 20

margin = 1
n_frames = 2
nx = 14
ny = 9

photons = [[10, 4000]]
pixel_size = 100.0

psf_size = 20
psf_z_range = 0.75

test_z_offset = 0.0
test_z_range = 0.100

tolerance = 0.3
x_size = 300
y_size = 200

z_planes = [0.1, -0.1, 0.1, -0.1]
z_value = [-0.3, 0.0, 0.3]
