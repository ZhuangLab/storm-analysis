#!/usr/bin/env python
"""
Settings to use in SLURM (multiplane) diagnostics. For this we just use the pupil function
PSF model as it is the simplest to setup.

Note: Background photons are per plane, total photons are divided across all the planes.

Hazen 08/18
"""
import numpy

camera_gain = 1.0
camera_offset = 100.0
camera_variance = 1.0
divisions = 10
independent_heights = 0
iterations = 20

# Mapping with a small offset.
if True:
    mappings = {"0_0_x" : numpy.array([0.0, 1.0, 0.0]),
                "0_0_y" : numpy.array([0.0, 0.0, 1.0]),
                "0_1_x" : numpy.array([2.0, 1.0, 0.0]),
                "0_1_y" : numpy.array([5.0, 0.0, 1.0]),
                "1_0_x" : numpy.array([-2.0, 1.0, 0.0]),
                "1_0_y" : numpy.array([-5.0, 0.0, 1.0])}

# Mapping with x flip.
if False:
    mappings = {"0_0_x" : numpy.array([0.0, 1.0, 0.0]),
                "0_0_y" : numpy.array([0.0, 0.0, 1.0]),
                "0_1_x" : numpy.array([302.0, -1.0, 0.0]),
                "0_1_y" : numpy.array([5.0, 0.0, 1.0]),
                "1_0_x" : numpy.array([302.0, -1.0, 0.0]),
                "1_0_y" : numpy.array([-5.0, 0.0, 1.0])}

margin = 1
n_frames = 100
nx = 14
ny = 9

off_time = 400.0
on_time = 2.0

photons = [10, 4000]
pixel_size = 100.0

psf_size = 30

pupil_fn = []

test_z_offset = 0.0
test_z_range = 0.0

tolerance = 0.3
x_size = 300
y_size = 200
z_planes = [-0.250, 0.250]
z_value = [-0.3, 0.0, 0.3]

pupilfn_z_range = 0.75

wdir = "slurm_test"
