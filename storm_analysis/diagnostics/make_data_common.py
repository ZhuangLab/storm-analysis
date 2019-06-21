#!/usr/bin/env python
"""
That which is common to all of the XYZ.make_data diagnostic modules.

Hazen 06/19
"""
import numpy
import os

import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate


def makeDataPupilFn(settings, dither):
    """
    Pupil function PSF, Ideal camera movies.
    """
    index = 1

    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, settings.camera_offset)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, settings.zmn)
        
        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                dither = dither,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
            
        sim.simulate(wdir + "/test.tif", "grid_list.hdf5", settings.n_frames)

        index += 1

    makePeakFile(settings)

    
def makeDataPupilFnCMOS(settings, cal_file, dither):
    """
    Pupil function PSF, CMOS camera movies.
    """
    index = 1

    for [bg, photons] in settings.photons:

        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)

        bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
        cam_f = lambda s, x, y, i3 : camera.SCMOS(s, x, y, i3, cal_file)
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, settings.zmn)
        
        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                dither = dither,
                                x_size = settings.x_size,
                                y_size = settings.y_size)
            
        sim.simulate(wdir + "/test.tif", "grid_list.hdf5", settings.n_frames)

        index += 1

    makePeakFile(settings)

    
def makePeakFile(settings):
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

