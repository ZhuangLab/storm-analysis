#!/usr/bin/env python
"""
Configure folder for Pupilfn testing.

Hazen 10/17
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
    Create a Pupilfn parameters object.
    """
    params = parameters.ParametersPupilFn()

    params.setAttr("max_frame", "int", -1) 
    params.setAttr("start_frame", "int", -1)
    
    params.setAttr("background_sigma", "float", 8.0)
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


# Create parameters file for analysis.
#
print("Creating XML file.")
params = testingParameters()
params.toXMLFile("pupilfn.xml")

# Create localization on a grid file.
#
print("Creating gridded localization.")
sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
subprocess.call(["python", sim_path + "emitters_on_grid.py",
                 "--bin", "grid_list.hdf5",
                 "--nx", str(settings.nx),
                 "--ny", str(settings.ny),
                 "--spacing", "20",
                 "--zrange", str(settings.test_z_range),
                 "--zoffset", str(settings.test_z_offset)])

# Create randomly located localizations file.
#
print("Creating random localization.")
subprocess.call(["python", sim_path + "emitters_uniform_random.py",
                 "--bin", "random_list.hdf5",
                 "--density", "1.0",
                 "--margin", str(settings.margin),
                 "--sx", str(settings.x_size),
                 "--sy", str(settings.y_size),
                 "--zrange", str(settings.test_z_range)])

# Create the pupil function.
#
pupilfn_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/pupilfn/"
print("Creating pupil function.")
subprocess.call(["python", pupilfn_path + "make_pupil_fn.py",
                 "--filename", "pupil_fn.pfn",
                 "--size", str(settings.pupil_size),
                 "--pixel-size", str(settings.pixel_size),
                 "--zmn", str(settings.zmn)])

                
