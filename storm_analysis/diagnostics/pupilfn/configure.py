#!/usr/bin/env python
"""
Configure folder for Pupilfn testing.

Hazen 06/19
"""
import numpy
import os

import storm_analysis
import storm_analysis.sa_library.parameters as parameters

import storm_analysis.simulator.emitters_on_grid as emittersOnGrid
import storm_analysis.simulator.emitters_uniform_random as emittersUniformRandom

import storm_analysis.pupilfn.make_pupil_fn as makePupilFn

import storm_analysis.diagnostics.pupilfn.settings as settings


def testingParameters(cal_file = None):
    """
    Create a Pupilfn parameters object.
    """
    params = parameters.ParametersPupilFn()

    params.setAttr("max_frame", "int", -1) 
    params.setAttr("start_frame", "int", -1)
    
    params.setAttr("background_sigma", "float", 8.0)

    if cal_file is not None:
        params.setAttr("camera_calibration", "filename", cal_file)
    else:
        params.setAttr("camera_gain", "float", settings.camera_gain)
        params.setAttr("camera_offset", "float", settings.camera_offset)
        
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("pupil_function", "filename", "pupil_fn.pfn")
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", 6.0)

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.0")

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    # Set Z range.
    params.setAttr("max_z", "float", 0.5)
    params.setAttr("min_z", "float", -0.5)

    # 'peak_locations' testing.
    if hasattr(settings, "peak_locations") and (settings.peak_locations is not None):
        params.setAttr("peak_locations", "filename", settings.peak_locations)

    return params


def configure(cal_file = None):
    
    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters(cal_file = cal_file)
    params.toXMLFile("pupilfn.xml")

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

    # Create the pupil function.
    #
    print("Creating pupil function.")
    makePupilFn.makePupilFunction("pupil_fn.pfn",
                                  settings.pupil_size,
                                  settings.pixel_size * 1.0e-3,
                                  settings.zmn)

                
if (__name__ == "__main__"):
    configure()
    
