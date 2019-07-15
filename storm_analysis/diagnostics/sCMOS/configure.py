#!/usr/bin/env python
"""
Configure folder for sCMOS testing.

Hazen 09/17
"""
import numpy
import os

import storm_analysis
import storm_analysis.sa_library.parameters as parameters

import storm_analysis.simulator.emitters_on_grid as emittersOnGrid
import storm_analysis.simulator.emitters_uniform_random as emittersUniformRandom

import storm_analysis.diagnostics.sCMOS.settings as settings


def testingParameters(cal_file):
    """
    Create a sCMOS parameters object.
    """
    params = parameters.ParametersSCMOS()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)
    params.setAttr("verbosity", "int", settings.verbosity)

    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("camera_calibration", "filename", cal_file)
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("fit_error_model", "string", settings.fit_error_model)
    params.setAttr("foreground_sigma", "float", 1.5)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("model", "string", settings.model)
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("roi_size", "int", settings.roi_size)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", settings.threshold)

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.0")

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    # Z fitting.
    #
    # These are nonsense values. We test either '2D' of '3D' mode
    # and check how well we do at fitting the localization widths.
    #
    params.setAttr("do_zfit", "int", 0)

    params.setAttr("cutoff", "float", 0.0)    
    params.setAttr("max_z", "float", 0.5)
    params.setAttr("min_z", "float", -0.5)
    params.setAttr("z_value", "float", 0.0)
    params.setAttr("z_step", "float", 1.0)

    params.setAttr("wx_wo", "float", 1.0)
    params.setAttr("wx_c", "float", 1.0)
    params.setAttr("wx_d", "float", 1.0)
    params.setAttr("wxA", "float", 0.0)
    params.setAttr("wxB", "float", 0.0)
    params.setAttr("wxC", "float", 0.0)
    params.setAttr("wxD", "float", 0.0)

    params.setAttr("wy_wo", "float", 1.0)
    params.setAttr("wy_c", "float", 1.0)
    params.setAttr("wy_d", "float", 1.0)
    params.setAttr("wyA", "float", 0.0)
    params.setAttr("wyB", "float", 0.0)
    params.setAttr("wyC", "float", 0.0)
    params.setAttr("wyD", "float", 0.0)

    # 'peak_locations' testing.
    if hasattr(settings, "peak_locations") and (settings.peak_locations is not None):
        params.setAttr("peak_locations", "filename", settings.peak_locations)

    return params

    
def configure(cal_file = None):

    # Create sCMOS calibration file if not specified.
    #
    if cal_file is None:
        cal_file = "calib.npy"
        offset = numpy.zeros((settings.y_size, settings.x_size)) + settings.camera_offset
        variance = numpy.ones((settings.y_size, settings.x_size)) * settings.camera_variance
        gain = numpy.ones((settings.y_size, settings.x_size)) * settings.camera_gain
        rqe = numpy.ones((settings.y_size, settings.x_size))
        numpy.save(cal_file, [offset, variance, gain, rqe, 2])

    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters(cal_file)
    params.toXMLFile("scmos.xml", pretty = True)

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
                                                10,
                                                settings.x_size,
                                                settings.y_size,
                                                0.0)

if (__name__ == "__main__"):
    configure()
    
