#!/usr/bin/env python
"""
Configure folder for CS decon testing.

Hazen 11/19
"""
import argparse
import numpy

import storm_analysis.sa_library.parameters as parameters

import storm_analysis.diagnostics.spliner_configure_common as splinerConfigureCommon

import storm_analysis.diagnostics.fista_decon.settings as settings


def testingParameters():
    """
    Create a Spliner parameters object for CS deconvolution.
    """
    if (settings.decon_method == "ADMM"):
        params = parameters.ParametersSpliner()
        params.update(parameters.ParametersADMM())
        params.update(parameters.ParametersRollingBall())
    elif (settings.decon_method == "3denseSTORM"):
        params = parameters.ParametersSpliner()
        params.update(parameters.Parameters3denseSTORM())
        params.update(parameters.ParametersRollingBall())
    elif (settings.decon_method == "FISTA"):
        params = parameters.ParametersSpliner()
        params.update(parameters.ParametersFISTA())
        params.update(parameters.ParametersRollingBall())
    else:
        params = parameters.ParametersSplinerSTD()
        
    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)    
    
    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("camera_gain", "float", settings.camera_gain)
    params.setAttr("camera_offset", "float", settings.camera_offset)
    params.setAttr("find_max_radius", "int", 15)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("max_z", "float", 1.0)
    params.setAttr("min_z", "float", -1.0)
    params.setAttr("no_fitting", "int", 1)

    params.setAttr("sigma", "float", 1.5)
    params.setAttr("spline", "filename", "psf.spline")
    params.setAttr("threshold", "float", 6.0)

    # ADMM.
    if (settings.decon_method == "ADMM"):
        params.setAttr("decon_method", "string", settings.decon_method)
        params.setAttr("admm_iterations", "int", settings.cs_iterations)
        params.setAttr("admm_lambda", "float", 10.0)
        params.setAttr("admm_number_z", "int", 5)
        params.setAttr("admm_rho", "float", 0.1)
        params.setAttr("admm_threshold", "float", 500.0)

        params.setAttr("background_estimator", "string", "RollingBall")
        params.setAttr("rb_radius", "float", 10.0)
        params.setAttr("rb_sigma", "float", 1.0)

    # 3denseSTORM
    elif (settings.decon_method == "3denseSTORM"):
        params.setAttr("decon_method", "string", settings.decon_method)
        params.setAttr("ds3_beta", "float", 2.44)
        params.setAttr("ds3_eta", "float", 9.0e-3)
        params.setAttr("ds3_iterations", "int", settings.cs_iterations)
        params.setAttr("ds3_micro", "float", 0.1)
        params.setAttr("ds3_number_z", "int", 5)
        params.setAttr("ds3_threshold", "float", 500.0)
        
        params.setAttr("background_estimator", "string", "RollingBall")
        params.setAttr("rb_radius", "float", 10.0)
        params.setAttr("rb_sigma", "float", 1.0)
        
    # FISTA
    elif (settings.decon_method == "FISTA"):
        params.setAttr("decon_method", "string", settings.decon_method)
        params.setAttr("fista_iterations", "int", settings.cs_iterations)
        params.setAttr("fista_lambda", "float", 20.0)
        params.setAttr("fista_number_z", "int", 5)
        params.setAttr("fista_threshold", "float", 500.0)
        params.setAttr("fista_timestep", "float", 0.25)
        
        params.setAttr("background_estimator", "string", "RollingBall")
        params.setAttr("rb_radius", "float", 10.0)
        params.setAttr("rb_sigma", "float", 1.0)

    # Standard
    else:
        params.setAttr("z_value", "float-array", [-0.6,-0.3,0.0,0.3,0.6])
        
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
    

def configure(no_splines):
    
    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters()
    params.toXMLFile("fdecon.xml")

    # Standard spliner configuration.
    #
    splinerConfigureCommon.configure(settings, settings.use_dh, no_splines)


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description = 'FISTA decon configuration.')
    
    parser.add_argument('--no-splines', dest='no_splines', action='store_true', default = False)
    
    args = parser.parse_args()
    
    configure(args.no_splines)


