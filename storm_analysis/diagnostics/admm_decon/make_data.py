#!/usr/bin/env python
"""
Make data for testing ADMM Spliner.

Hazen 02/18
"""
import numpy
import os

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.diagnostics.make_data_common as makeDataCommon

import storm_analysis.diagnostics.admm_decon.settings as settings


def makeData(dither = False):
    """
    Ideal camera movies.
    """
    assert not settings.use_dh    
    makeDataCommon.makeDataPupilFn(settings, dither)

    
def makeDataDoubleHelix(dither = False):
    """
    Make data simulating Double Helix PSF.
    """
    assert settings.use_dh

    index = 1

    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.DHPSF(s, x, y, i3, 100.0, z_range = settings.spline_z_range)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
    
        sim.simulate(wdir + "/test.dax", "grid_list.hdf5", settings.n_frames)
        
        index += 1

        
if (__name__ == "__main__"):
    makeData()
