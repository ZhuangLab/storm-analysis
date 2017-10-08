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
# Total analysis time 331.59 seconds
# Recall 0.52341
# Noise 0.46885
# XYZ Error (nm):
# test_01	35.05	35.32	74.93
# test_02	18.62	18.25	43.44
#
# independent_heights = 1
#
# Analysis Summary:
# Total analysis time 308.70 seconds
# Recall 0.53548
# Noise 0.46109
# XYZ Error (nm):
# test_01	34.20	34.41	69.40
# test_02	18.33	17.96	40.87
#

if True:

    # Create .bin files for each plane.
    i3_locs = readinsight3.loadI3File("grid_list.bin")

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
        pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, photons)
        psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, [])

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

