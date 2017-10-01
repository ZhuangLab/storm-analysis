#!/usr/bin/env python
"""
Make data for testing sCMOS. 

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

# sCMOS calibration file.
numpy.save("calib.npy", [numpy.zeros((xsize, ysize)) + 100.0,  # offset.
                         numpy.ones((xsize, ysize))*16.0,      # variance (ADU^2).
                         numpy.ones((xsize, ysize))*2.0])      # gain (ADU/e-)

# sCMOS camera movies.
#
# For these simulations we expect (approximately) these results:
#
# Analysis Summary:
# Total analysis time 10.64 seconds
# Recall 0.93726
# Noise 0.05972
# XY Error (nm):
# test_01	14.44	14.53
# test_02	8.26	8.24
#
# XY Width Error, Mean difference with truth, Standard deviation (pixels):
# test_01	0.029	0.116	0.029	0.116
# test_02	0.017	0.105	0.017	0.105
#
if True:
    for [bg, photons] in [[20, 500], [20, 1000]]:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        sim = simulate.Simulate(background_factory = lambda settings, xs, ys, i3data : background.UniformBackground(settings, xs, ys, i3data, photons = bg),
                                camera_factory = lambda settings, xs, ys, i3data : camera.SCMOS(settings, xs, ys, i3data, 0.0, "calib.npy"),
                                photophysics_factory = lambda settings, xs, ys, i3data : photophysics.AlwaysOn(settings, xs, ys, i3data, photons = photons),
                                psf_factory = lambda settings, xs, ys, i3data : psf.GaussianPSF(settings, xs, ys, i3data, pixel_size),
                                x_size = xsize, y_size = ysize)
    
        sim.simulate(wdir + "/test.dax", "grid_list.bin", n_frames)
        
        index += 1
