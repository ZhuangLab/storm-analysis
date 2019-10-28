#!/usr/bin/env python
"""
Configure folder for Multiplane 3D-DAOSTORM testing.

Hazen 01/18
"""
import argparse
import numpy
import os
import pickle

import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.emitters_on_grid as emittersOnGrid
import storm_analysis.simulator.emitters_uniform_random as emittersUniformRandom

import storm_analysis.diagnostics.multiplane_dao.settings as settings


def testingParameters():
    """
    Create a Multiplane 3D-DAOSTORM parameters object.
    """
    params = parameters.ParametersMultiplaneDao()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)
    
    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("find_max_radius", "int", 2)
    params.setAttr("foreground_sigma", "float", 1.5)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("mapping", "filename", "map.map")
    params.setAttr("no_fitting", "int", 0)
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("roi_size", "int", settings.roi_size)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", 6.0)

    params.setAttr("channel0_cal", "filename", "calib.npy")
    params.setAttr("channel1_cal", "filename", "calib.npy")

    params.setAttr("channel0_ext", "string", "_c1.dax")
    params.setAttr("channel1_ext", "string", "_c2.dax")

    params.setAttr("channel0_offset", "int", 0)
    params.setAttr("channel1_offset", "int", 0)

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.0")
    
    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    # 'peak_locations' testing.
    if hasattr(settings, "peak_locations") and (settings.peak_locations is not None):
        params.setAttr("peak_locations", "filename", settings.peak_locations)    

    return params
    
def configure():
    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters()
    params.toXMLFile("multiplane.xml")

    # Create localization on a grid file.
    #
    print("Creating gridded localization.")
    emittersOnGrid.emittersOnGrid("grid_list.hdf5",
                                  settings.nx,
                                  settings.ny,
                                  1.5,
                                  20,
                                  settings.test_z_range,
                                  settings.test_z_offset)

    # Create randomly located localizations file.
    #
    print("Creating random localization.")
    emittersUniformRandom.emittersUniformRandom("random_list.hdf5",
                                                1.0,
                                                settings.margin,
                                                settings.x_size,
                                                settings.y_size,
                                                settings.test_z_range)
    
    # Create sCMOS camera calibration files.
    #
    numpy.save("calib.npy", [numpy.zeros((settings.y_size, settings.x_size)) + settings.camera_offset,
                             numpy.ones((settings.y_size, settings.x_size)) * settings.camera_variance,
                             numpy.ones((settings.y_size, settings.x_size)) * settings.camera_gain,
                             numpy.ones((settings.y_size, settings.x_size)),
                             2])

    # Create mapping file.
    with open("map.map", 'wb') as fp:
        pickle.dump(settings.mappings, fp)


if (__name__ == "__main__"):
    configure()
    
