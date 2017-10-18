#!/usr/bin/env python
"""
Configure folder for PSF FFT testing.

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
    Create a PSF FFT parameters object.
    """
    params = parameters.ParametersPSFFFT()

    params.setAttr("max_frame", "int", -1) 
    params.setAttr("start_frame", "int", -1)
    params.setAttr("append_metadata", "int", 0)
    
    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("camera_gain", "float", settings.camera_gain)
    params.setAttr("camera_offset", "float", settings.camera_offset)
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("iterations", "int", 20)
    params.setAttr("orientation", "string", "normal")
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("psf", "filename", "psf.psf")
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

    # Use pre-specified fitting locations.
    if False:
        params.setAttr("peak_locations", "filename", "olist.bin")

    return params


# Create parameters file for analysis.
#
print("Creating XML file.")
params = testingParameters()
params.toXMLFile("psf_fft.xml")

# Create localization on a grid file.
#
print("Creating gridded localization.")
sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
subprocess.call(["python", sim_path + "emitters_on_grid.py",
                 "--bin", "grid_list.bin",
                 "--nx", "14",
                 "--ny", "9",
#                 "--nx", "1",
#                 "--ny", "1",
                 "--spacing", "20",
                 "--zrange", str(settings.test_z_range)])

# Create randomly located localizations file.
#
print("Creating random localization.")
subprocess.call(["python", sim_path + "emitters_uniform_random.py",
                 "--bin", "random_list.bin",
                 "--density", "1.0",
                 "--sx", str(settings.x_size),
                 "--sy", str(settings.y_size)])

# Create the PSF.
#
pupilfn_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/psf_fft/"
print("Creating (theoritical) psf.")
subprocess.call(["python", pupilfn_path + "make_psf_from_pf.py",
                 "--filename", "psf.psf",
                 "--size", str(settings.psf_size),
                 "--pixel-size", str(settings.psf_size),
                 "--zrange", str(settings.psf_z_range)])
