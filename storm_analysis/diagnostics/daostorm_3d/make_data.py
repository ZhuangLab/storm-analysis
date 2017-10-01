#!/usr/bin/env python
"""
Make data for testing 3D-DAOSTORM. 

The default tests are pretty easy as they are just relatively bright
localizations on a grid.

Hazen 09/17
"""
import numpy
import os

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import settings

index = 1

# Ideal camera movies.
#
# For these simulations we expect (approximately) these results:
#
# Analysis Summary:
# Total analysis time 9.13 seconds
# Recall 0.95278
# Noise 0.04597
# XY Error (nm):
# test_01	13.78	13.80
# test_02	7.84	7.86
#
# XY Width Error, Mean difference with truth, Standard deviation (pixels):
# test_01	0.034	0.112	0.034	0.112
# test_02	0.031	0.098	0.031	0.098
#
if True:
    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, [[1.3, 2, 2]])
        #psf_f = lambda s, xs, ys, i3 : psf.GaussianPSF(settings, xs, ys, i3data, settings.pixel_size),

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
    
        sim.simulate(wdir + "/test.dax", "grid_list.bin", settings.n_frames)
        
        index += 1
        
