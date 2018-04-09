#!/usr/bin/env python
"""
Make data for testing Multiplane. 

The default tests are pretty easy as they are just relatively bright
localizations on a grid.

Hazen 09/17
"""
import numpy
import os
import pickle

import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.diagnostics.multicolor.settings as settings


def makeData():
    index = 1

    # Create HDF5 files for each plane.
    #
    for elt in ["grid_list.hdf5", "random_storm.hdf5"]:
        locs = saH5Py.loadLocalizations(elt)
        locs["color"] = numpy.random.randint(4, size = locs["x"].size)
        zo = locs["z"].copy()
        
        locs["z"][:] = zo + 1.0e-3 * settings.z_planes[0]
        saH5Py.saveLocalizations("sim_input_c1_" + elt, locs)
        for i in range(1,4):
            locs["x"] += settings.dx
            locs["y"] += settings.dy
            locs["z"][:] = zo + 1.0e-3 * settings.z_planes[i]
            saH5Py.saveLocalizations("sim_input_c" + str(i+1) + "_" + elt, locs)

    if True:
        
        # Create a movie for each plane.
        for [bg, photons] in settings.photons:

            # Adjust photons by the number of planes.
            photons = photons/float(len(settings.z_planes))

            wdir = "test_{0:02d}".format(index)
            print(wdir)
            if not os.path.exists(wdir):
                os.makedirs(wdir)

            for i in range(4):
                bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
                cam_f = lambda s, x, y, i3 : camera.SCMOS(s, x, y, i3, "calib.npy")
                pp_f = lambda s, x, y, i3 : photophysics.AlwaysOnMC(s, x, y, i3, color = i, photons = photons)
                psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, [])

                sim = simulate.Simulate(background_factory = bg_f,
                                        camera_factory = cam_f,
                                        photophysics_factory = pp_f,
                                        psf_factory = psf_f,
                                        x_size = settings.x_size,
                                        y_size = settings.y_size)

            
                sim.simulate(wdir + "/test_c" + str(i+1) + ".dax",
                             "sim_input_c" + str(i+1) + "_grid_list.hdf5",
                             settings.n_frames)
        
            index += 1


if (__name__ == "__main__"):
    makeData()
    


