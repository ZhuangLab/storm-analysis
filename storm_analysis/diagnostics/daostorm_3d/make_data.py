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
# Total analysis time 9.60 seconds
# Recall 0.95139
# Noise 0.04702
# XY Error (nm):
# test_01	13.60	13.87
# test_02	7.89	7.87
# 
# XY Width Error, Mean difference with truth, Standard deviation (pixels):
# test_01	0.028	0.106	0.028	0.106
# test_02	0.019	0.101	0.019	0.101
#

# Gaussian PSF, uniform background.
if True:
    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.GaussianPSF(s, x, y, i3, settings.pixel_size)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
    
        sim.simulate(wdir + "/test.dax", "grid_list.bin", settings.n_frames)
        
        index += 1

# Pupil Function PSF.
if False:
    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, [])

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
    
        sim.simulate(wdir + "/test.dax", "grid_list.bin", settings.n_frames)
        
        index += 1        

# Gaussian non-uniform background, always on.
if False:
    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.GaussianBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.GaussianPSF(s, x, y, i3, settings.pixel_size)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
    
        sim.simulate(wdir + "/test.dax", "grid_list.bin", settings.n_frames)
        
        index += 1

# Gaussian non-uniform background, STORM.
if False:
    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.GaussianBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
        pp_f = lambda s, x, y, i3 : photophysics.SimpleSTORM(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.GaussianPSF(s, x, y, i3, settings.pixel_size)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
    
        sim.simulate(wdir + "/test.dax", "random_list.bin", settings.n_frames)
        
        index += 1        

# Sloped non-uniform background, always on.
if False:
    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.SlopedBackground(s, x, y, i3, slope = 0.4, offset = 10)
        cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.GaussianPSF(s, x, y, i3, settings.pixel_size)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                dither = True,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
    
        sim.simulate(wdir + "/test.dax", "grid_list.bin", settings.n_frames)
        
        index += 1

# Sine non-uniform background, always on.
if False:
    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.SineBackground(s, x, y, i3, photons = bg, period = 45)
        cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.GaussianPSF(s, x, y, i3, settings.pixel_size)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                dither = True,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
    
        sim.simulate(wdir + "/test.dax", "grid_list.bin", settings.n_frames)
        
        index += 1
        
