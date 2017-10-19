#!/usr/bin/env python
"""
Configure folder for Spliner testing.

Hazen 09/17
"""
import argparse
import inspect
import numpy
import os
import subprocess

import storm_analysis
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.readinsight3 as readinsight3

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import settings


parser = argparse.ArgumentParser(description = 'Spline diagnostics configuration.')

parser.add_argument('--no-splines', dest='no_splines', action='store_true', default = False)

args = parser.parse_args()
    

def testingParameters():
    """
    Create a Spliner parameters object.
    """
    params = parameters.ParametersSpliner()

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

    # Use pre-specified fitting locations.
    if False:
        params.setAttr("peak_locations", "filename", "olist.bin")
        
    return params
    

# Create parameters file for analysis.
#
print("Creating XML file.")
params = testingParameters()
params.toXMLFile("spliner.xml")

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
                 "--sy", str(settings.y_size),
                 "--zrange", str(settings.test_z_range)])

# Create sparser grid for PSF measurement.
#
print("Creating data for PSF measurement.")
sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
subprocess.call(["python", sim_path + "emitters_on_grid.py",
                 "--bin", "sparse_list.bin",
                 "--nx", "6",
                 "--ny", "3",
                 "--spacing", "40"])


if args.no_splines:
    exit()

    
# Create beads.txt file for spline measurement.
#
locs = readinsight3.loadI3File("sparse_list.bin")
numpy.savetxt("beads.txt", numpy.transpose(numpy.vstack((locs['xc'], locs['yc']))))

# Create drift file, this is used to displace the localizations in the
# PSF measurement movie.
#
dz = numpy.arange(-settings.spline_z_range, settings.spline_z_range + 5.0, 10.0)
drift_data = numpy.zeros((dz.size, 3))
drift_data[:,2] = dz
numpy.savetxt("drift.txt", drift_data)

# Also create the z-offset file.
#
z_offset = numpy.ones((dz.size, 2))
z_offset[:,1] = dz
numpy.savetxt("z_offset.txt", z_offset)

# Create simulated data for PSF measurement.
#
bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = 10)
cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, 100.)
drift_f = lambda s, x, y, i3 : drift.DriftFromFile(s, x, y, i3, "drift.txt")
pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, 20000.0)
psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, 100.0, settings.zmn)

sim = simulate.Simulate(background_factory = bg_f,
                        camera_factory = cam_f,
                        drift_factory = drift_f,
                        photophysics_factory = pp_f,
                        psf_factory = psf_f,
                        x_size = settings.x_size,
                        y_size = settings.y_size)
                        
sim.simulate("spline.dax", "sparse_list.bin", dz.size)

# Measure the PSF.
#
print("Measuring PSF.")
spliner_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/spliner/"
subprocess.call(["python", spliner_path + "measure_psf_beads.py",
                 "--movie", "spline.dax",
                 "--zoffset", "z_offset.txt",
                 "--beads", "beads.txt",
                 "--psf", "psf.psf"])

# Measure the Spline.
#
print("Measuring Spline.")
subprocess.call(["python", spliner_path + "psf_to_spline.py",
                 "--psf", "psf.psf",
                 "--spline", "psf.spline",
                 "--spline_size", str(settings.spline_size)])
