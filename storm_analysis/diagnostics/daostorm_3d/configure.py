#!/usr/bin/env python
"""
Configure folder for 3D-DAOSTORM testing.

Hazen 09/17
"""
import inspect
import numpy
import os

import storm_analysis
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.emitters_on_grid as emittersOnGrid
import storm_analysis.simulator.emitters_uniform_random as emittersUniformRandom

import storm_analysis.diagnostics.daostorm_3d.settings as settings


def testingParameters():
    """
    Create a 3D-DAOSTORM parameters object.
    """
    params = parameters.ParametersDAO()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)

    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("camera_gain", "float", settings.camera_gain)
    params.setAttr("camera_offset", "float", settings.camera_offset)
    params.setAttr("fftw_estimate", "int", settings.fftw_estimate)
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("fit_error_model", "string", settings.fit_error_model)
    params.setAttr("foreground_sigma", "float", 1.0)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("model", "string", settings.model)
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("roi_size", "int", 9)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", settings.threshold)

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0")

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    # Z fitting.
    #
    # These are nonsense values. We test either '2D' or '3D' mode
    # and check how well we do at fitting the localization widths.
    #
    params.setAttr("do_zfit", "int", 0)

    params.setAttr("cutoff", "float", 0.0)    
    params.setAttr("max_z", "float", 0.5)
    params.setAttr("min_z", "float", -0.5)
    params.setAttr("z_value", "float", 0.0)
    params.setAttr("z_step", "float", 1.0)

    params.setAttr("wx_wo", "float", 300.0)
    params.setAttr("wx_c", "float", 150.0)
    params.setAttr("wx_d", "float", 400.0)
    params.setAttr("wxA", "float", 0.0)
    params.setAttr("wxB", "float", 0.0)
    params.setAttr("wxC", "float", 0.0)
    params.setAttr("wxD", "float", 0.0)

    params.setAttr("wy_wo", "float", 300.0)
    params.setAttr("wy_c", "float", -150.0)
    params.setAttr("wy_d", "float", 400.0)
    params.setAttr("wyA", "float", 0.0)
    params.setAttr("wyB", "float", 0.0)
    params.setAttr("wyC", "float", 0.0)
    params.setAttr("wyD", "float", 0.0)

    # File conversion testing.
    #params.setAttr("convert_to", "string", ".bin,.txt")
    #params.setAttr("convert_to", "string", ".bin")
    #params.setAttr("convert_to", "string", ".txt")

    # peak initialization testing.
    if hasattr(settings, "no_fitting") and settings.no_fitting:
        params.setAttr("no_fitting", "int", 1)
        
    # 'peak_locations' testing.
    if hasattr(settings, "peak_locations") and (settings.peak_locations is not None):
        params.setAttr("peak_locations", "filename", settings.peak_locations)

    # mask testing.
    if hasattr(settings, "x_start") and (settings.x_start is not None):
        params.setAttr("x_start", "int", settings.x_start)

    if hasattr(settings, "x_stop") and (settings.x_stop is not None):
        params.setAttr("x_stop", "int", settings.x_stop)

    if hasattr(settings, "y_start") and (settings.y_start is not None):
        params.setAttr("y_start", "int", settings.y_start)

    if hasattr(settings, "y_stop") and (settings.y_stop is not None):
        params.setAttr("y_stop", "int", settings.y_stop)        

    return params
    
def configure():
    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters()
    params.toXMLFile("dao.xml")

    # Create localization on a grid file.
    #
    print("Creating gridded localization.")
    emittersOnGrid.emittersOnGrid("grid_list.hdf5",
                                  settings.nx,
                                  settings.ny,
                                  1.5,
                                  20,
                                  0.0,
                                  0.0)
    
    # Create randomly located localizations file.
    #
    print("Creating random localization.")
    emittersUniformRandom.emittersUniformRandom("random_list.hdf5",
                                                1.0,
                                                settings.margin,
                                                settings.x_size,
                                                settings.y_size,
                                                0.0)


if (__name__ == "__main__"):
    configure()
    
