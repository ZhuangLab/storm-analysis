#!/usr/bin/env python
"""
Make data for testing drift correction.

Hazen 01/18
"""
import numpy
import os

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.diagnostics.drift_correction.settings as settings


def makeData():
    index = 1
    
    # Gaussian non-uniform background, STORM.
    if True:

        for elt in ["drift_xy.txt", "drift_xyz.txt"]:
            bg = settings.background
            photons = settings.photons
    
            wdir = "test_{0:02d}".format(index)
            print(wdir)
            if not os.path.exists(wdir):
                os.makedirs(wdir)

            bg_f = lambda s, x, y, i3 : background.GaussianBackground(s, x, y, i3, photons = bg)
            cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
            drift_f = lambda s, x, y, i3 : drift.DriftFromFile(s, x, y, i3, elt)
            pp_f = lambda s, x, y, i3 : photophysics.SimpleSTORM(s, x, y, i3, photons)
            psf_f = lambda s, x, y, i3 : psf.GaussianPSF(s, x, y, i3, settings.pixel_size)
        
            sim = simulate.Simulate(background_factory = bg_f,
                                    camera_factory = cam_f,
                                    drift_factory = drift_f,
                                    photophysics_factory = pp_f,
                                    psf_factory = psf_f,
                                    x_size = settings.x_size,
                                    y_size = settings.y_size)
        
            sim.simulate(wdir + "/test.dax", "clusters_list.hdf5", settings.n_frames)
            #sim.simulate(wdir + "/test.dax", "lines_list.hdf5", settings.n_frames)
        
            index += 1


if (__name__ == "__main__"):
    makeData()
    

