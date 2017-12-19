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

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import settings

index = 1

# sCMOS camera movies.
#
# For these simulations we expect (approximately) these results:
#
# independent_heights = 0
#
# Analysis Summary:
# Processed 24797 localizations in 307.92 seconds, 80.53/sec
# Recall 0.59452
# Noise 0.39581
# XYZ Error (nm):
# test_01	31.07	31.29	52.26
# test_02	16.63	16.62	28.57
#
# independent_heights = 1
#
# Analysis Summary:
# Processed 24825 localizations in 276.57 seconds, 89.76/sec
# Recall 0.59313
# Noise 0.39791
# XYZ Error (nm):
# test_01	31.49	31.71	80.12
# test_02	16.71	16.74	47.36
#

if True:

    grid = False

    # Create .bin files for each plane.
    if grid:
        i3_locs = readinsight3.loadI3File("grid_list.bin")
    else:
        i3_locs = readinsight3.loadI3File("random_list.bin")
        
    # Load channel to channel mapping file.
    with open("map.map", 'rb') as fp:
        mappings = pickle.load(fp)

    for i, z_plane in enumerate(settings.z_planes):
        cx = mappings["0_" + str(i) + "_x"]
        cy = mappings["0_" + str(i) + "_y"]
        i3_temp = i3_locs.copy()
        xi = i3_temp["x"]
        yi = i3_temp["y"]
        zi = i3_temp["z"]
        xf = cx[0] + cx[1] * xi + cx[2] * yi
        yf = cy[0] + cy[1] * xi + cy[2] * yi
        zf = zi + z_plane
        i3dtype.posSet(i3_temp, "x", xf)
        i3dtype.posSet(i3_temp, "y", yf)
        i3dtype.posSet(i3_temp, "z", zf)
        with writeinsight3.I3Writer("sim_input_c" + str(i+1) + ".bin") as i3w:
            i3w.addMolecules(i3_temp)
        
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

        if grid:
            pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        else:
            pp_f = lambda s, x, y, i3 : photophysics.SimpleSTORM(s, x, y, i3, photons)

        psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, settings.pupil_fn)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)

        for i in range(len(settings.z_planes)):
            sim.simulate(wdir + "/test_c" + str(i+1) + ".dax",
                         "sim_input_c" + str(i+1) + ".bin",
                         settings.n_frames)
        
        index += 1

