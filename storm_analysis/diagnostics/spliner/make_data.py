#!/usr/bin/env python
"""
Make data for testing Spliner. 

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

import storm_analysis.diagnostics.spliner.settings as settings


def makeData(dither = False):
    """
    Ideal camera movies.
    """
    makeDataCommon.makeDataPupilFn(settings, dither)


def makeDataRQE(cal_file = "calib.npy", dither = False):
    """
    RQE testing camera movies.
    """
    makeDataCommon.makeDataPupilFnCMOS(settings, cal_file, dither)

    
def makeDataSpliner():

    index = 1

    # Non-uniform background, STORM.
    #
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
            
            sim.simulate(wdir + "/test.dax", "random_list.hdf5", settings.n_frames)
        
            index += 1

    # Ideal camera movies, PSF using the measured spline.
    #
    if False:
        for [bg, photons] in settings.photons:
            
            wdir = "test_{0:02d}".format(index)
            print(wdir)
            if not os.path.exists(wdir):
                os.makedirs(wdir)
            
            bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
            cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
            pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
            psf_f = lambda s, x, y, i3 : psf.Spline(s, x, y, i3, settings.pixel_size, "psf.spline")

            sim = simulate.Simulate(background_factory = bg_f,
                                    camera_factory = cam_f,
                                    photophysics_factory = pp_f,
                                    psf_factory = psf_f,
                                    x_size = settings.x_size,
                                    y_size = settings.y_size)
    
            sim.simulate(wdir + "/test.dax", "grid_list.hdf5", settings.n_frames)
        
            index += 1
        
    # Create "peak_locations" file if needed.
    #
    if hasattr(settings, "peak_locations") and (settings.peak_locations is not None):
        with saH5Py.SAH5Py("test_01/test_ref.hdf5") as h5:
            locs = h5.getLocalizationsInFrame(0)
            
        if settings.peak_locations.endswith(".hdf5"):
            saH5Py.saveLocalizations(settings.peak_locations, locs)
        else:
            numpy.savetxt(settings.peak_locations,
                          numpy.transpose(numpy.vstack((locs['x'],
                                                        locs['y'],
                                                        locs['height'],
                                                        locs['background']))))

if (__name__ == "__main__"):
    makeData()
    
