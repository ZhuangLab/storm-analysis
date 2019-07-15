#!/usr/bin/env python
"""
Make data for testing sCMOS. 

The default tests are pretty easy as they are just relatively bright
localizations on a grid.

Hazen 09/17
"""
import numpy
import os

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.diagnostics.make_data_common as makeDataCommon

import storm_analysis.diagnostics.sCMOS.settings as settings


def makeData(cal_file = "calib.npy", dither = False, verbosity = 1):
    index = 1

    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.SCMOS(s, x, y, i3, cal_file)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.GaussianPSF(s, x, y, i3, settings.pixel_size)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                dither = dither,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
            
        sim.simulate(wdir + "/test.tif", "grid_list.hdf5", settings.n_frames, verbosity = verbosity)
        
        index += 1

    makeDataCommon.makePeakFile(settings)
    

if (__name__ == "__main__"):
    makeData()
    
