#!/usr/bin/env python
"""
Make data for testing Multiplane 3D-DAOSTORM 

The default tests are pretty easy as they are just relatively bright
localizations on a grid.

Hazen 01/18
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

import storm_analysis.diagnostics.multiplane_dao.settings as settings


def makeData():
    index = 1

    if True:

        # Create localization files for each plane.
        h5_locs = saH5Py.loadLocalizations("grid_list.hdf5")
        
        # Load channel to channel mapping file.
        with open("map.map", 'rb') as fp:
            mappings = pickle.load(fp)

        for i, z_plane in enumerate(settings.z_planes):
            cx = mappings["0_" + str(i) + "_x"]
            cy = mappings["0_" + str(i) + "_y"]
            xi = h5_locs["x"].copy()
            yi = h5_locs["y"].copy()
            zi = h5_locs["z"].copy()
            xf = cx[0] + cx[1] * xi + cx[2] * yi
            yf = cy[0] + cy[1] * xi + cy[2] * yi
            zf = zi + z_plane
            h5_temp = {"x" : xf,
                       "y" : yf,
                       "z" : zf}
            saH5Py.saveLocalizations("sim_input_c" + str(i+1) + ".hdf5", h5_temp)
        
        # Create a movie for each plane.
        for [bg, photons] in settings.photons:

            # Adjust photons by the number of planes.
            photons = photons/float(len(settings.z_planes))

            wdir = "test_{0:02d}".format(index)
            print(wdir)
            if not os.path.exists(wdir):
                os.makedirs(wdir)

            bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
            cam_f = lambda s, x, y, i3 : camera.SCMOS(s, x, y, i3, "calib.npy")
            pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
            psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, settings.pupil_fn)
            
            sim = simulate.Simulate(background_factory = bg_f,
                                    camera_factory = cam_f,
                                    photophysics_factory = pp_f,
                                    psf_factory = psf_f,
                                    x_size = settings.x_size,
                                    y_size = settings.y_size)

            for i in range(len(settings.z_planes)):
                sim.simulate(wdir + "/test_c" + str(i+1) + ".dax",
                             "sim_input_c" + str(i+1) + ".hdf5",
                             settings.n_frames)
        
            index += 1


    # Create "peak_locations" file if needed.
    #
    if hasattr(settings, "peak_locations") and (settings.peak_locations is not None):
        with saH5Py.SAH5Py("test_01/test_c1_ref.hdf5") as h5:
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
    
