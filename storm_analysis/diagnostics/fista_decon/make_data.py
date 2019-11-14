#!/usr/bin/env python
"""
Make data for testing FISTA Spliner.

The default tests are pretty easy as they are just relatively bright
localizations on a grid.

Hazen 11/17
"""
import numpy
import os

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.diagnostics.fista_decon.settings as settings


def makeData():

    index = 1
    
    # Ideal camera movies.
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
            if settings.use_dh:
                psf_f = lambda s, x, y, i3 : psf.DHPSF(s, x, y, i3, 100.0, z_range = settings.spline_z_range)
            else:
                psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, 100.0, settings.zmn)

            sim = simulate.Simulate(background_factory = bg_f,
                                    camera_factory = cam_f,
                                    photophysics_factory = pp_f,
                                    psf_factory = psf_f,
                                    x_size = settings.x_size,
                                    y_size = settings.y_size)
    
            sim.simulate(wdir + "/test.tif", "grid_list.hdf5", settings.n_frames)

            index += 1


if (__name__ == "__main__"):
    makeData()
    

