#!/usr/bin/env python
"""
Configure folder for 3D-DAOSTORM testing.

Hazen 09/17
"""
import inspect
import numpy
import os
import subprocess

import storm_analysis
import storm_analysis.sa_library.parameters as parameters

import settings

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
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("foreground_sigma", "float", 1.0)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("model", "string", settings.model)
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("roi_size", "int", 9)
    params.setAttr("sigma", "float", 150.0/settings.pixel_size)
    params.setAttr("threshold", "float", 6.0)

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.1")

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    # Z fitting.
    #
    # These are nonsense values. We test either '2D' of '3D' mode
    # and check how will we do at fitting the localization widths.
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
    
    return params
    

# Create parameters file for analysis.
#
print("Creating XML file.")
params = testingParameters()
params.toXMLFile("dao.xml")

# Create localization on a grid file.
#
print("Creating gridded localization.")
sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
subprocess.call(["python", sim_path + "emitters_on_grid.py",
                 "--bin", "grid_list.hdf5",
                 "--nx", str(settings.nx),
                 "--ny", str(settings.ny),
                 "--spacing", "20"])

# Create randomly located localizations file.
#
print("Creating random localization.")
subprocess.call(["python", sim_path + "emitters_uniform_random.py",
                 "--bin", "random_list.hdf5",
                 "--density", "1.0",
                 "--margin", str(settings.margin),
                 "--sx", str(settings.x_size),
                 "--sy", str(settings.y_size)])
