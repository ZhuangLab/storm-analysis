#!/usr/bin/env python
"""
Configure folder for Spliner testing.

Hazen 09/17
"""
import argparse
import numpy

import storm_analysis.sa_library.parameters as parameters

import storm_analysis.diagnostics.spliner_configure_common as splinerConfigureCommon

import storm_analysis.diagnostics.spliner.settings as settings


def testingParameters(cal_file = None):
    """
    Create a Spliner parameters object.
    """
    params = parameters.ParametersSpliner()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)    

    params.setAttr("background_sigma", "float", 8.0)

    if cal_file is not None:
        params.setAttr("camera_calibration", "filename", cal_file)
    else:
        params.setAttr("camera_gain", "float", settings.camera_gain)
        params.setAttr("camera_offset", "float", settings.camera_offset)

    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("fit_error_model", "string", settings.fit_error_model)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("no_fitting", "int", 0)
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("spline", "filename", "psf.spline")
    params.setAttr("threshold", "float", 6.0)

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


def configure(no_splines, cal_file = None):

    # Create sCMOS calibration file if requested.
    #
    if cal_file is not None:
        offset = numpy.zeros((settings.y_size, settings.x_size)) + settings.camera_offset
        variance = numpy.ones((settings.y_size, settings.x_size))
        gain = numpy.ones((settings.y_size, settings.x_size)) * settings.camera_gain
        rqe = numpy.ones((settings.y_size, settings.x_size))
        numpy.save(cal_file, [offset, variance, gain, rqe, 2])
        
    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters(cal_file = cal_file)
    params.toXMLFile("spliner.xml")

    # Standard spliner configuration.
    #
    splinerConfigureCommon.configure(settings, False, no_splines)
    

if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description = 'Spline diagnostics configuration.')

    parser.add_argument('--no-splines', dest='no_splines', action='store_true', default = False)

    args = parser.parse_args()
    
    configure(args.no_splines)
