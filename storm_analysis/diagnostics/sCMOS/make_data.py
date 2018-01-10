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

import settings

index = 1

# sCMOS calibration file.
if settings.random_variance:
    variance = numpy.random.exponential(scale = settings.camera_variance, size = (settings.y_size, settings.x_size))
    offset = settings.camera_offset + variance
    numpy.save("calib.npy", [offset,
                             variance,
                             numpy.ones((settings.y_size, settings.x_size)) * settings.camera_gain,
                             1])
else:
    numpy.save("calib.npy", [numpy.zeros((settings.y_size, settings.x_size)) + settings.camera_offset,
                             numpy.ones((settings.y_size, settings.x_size)) * settings.camera_variance,
                             numpy.ones((settings.y_size, settings.x_size)) * settings.camera_gain,
                             1])

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
    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.SCMOS(s, x, y, i3, "calib.npy")
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.GaussianPSF(s, x, y, i3, settings.pixel_size)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
    
        sim.simulate(wdir + "/test.tif", "grid_list.hdf5", settings.n_frames)
        
        index += 1
