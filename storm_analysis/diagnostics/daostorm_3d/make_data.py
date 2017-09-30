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

pixel_size = 100.0
n_frames = 100
xsize = 300
ysize = 200
index = 1

# Ideal camera movies.
#
if True:
    for [bg, photons] in [[20, 500], [20, 1000]]:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        sim = simulate.Simulate(background_factory = lambda settings, xs, ys, i3data : background.UniformBackground(settings, xs, ys, i3data, photons = bg),
                                camera_factory = lambda settings, xs, ys, i3data : camera.Ideal(settings, xs, ys, i3data, 100.0),
                                photophysics_factory = lambda settings, xs, ys, i3data : photophysics.AlwaysOn(settings, xs, ys, i3data, photons = photons),
                                psf_factory = lambda settings, xs, ys, i3data : psf.GaussianPSF(settings, xs, ys, i3data, pixel_size),
                                x_size = xsize, y_size = ysize)
    
        sim.simulate(wdir + "/test.dax", "grid_list.bin", n_frames)
        
        index += 1

