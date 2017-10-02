#!/usr/bin/env python
"""
Configure folder for Multiplane testing.

Hazen 10/17
"""
import inspect
import numpy
import os
import pickle
import subprocess

import storm_analysis
import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import settings


def testingParameters():
    """
    Create a Spliner parameters object.
    """
    params = parameters.ParametersMultiplane()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)    
    params.setAttr("append_metadata", "int", 0)
    
    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("find_max_radius", "int", 2)
    params.setAttr("independent_heights", "int", settings.independent_heights)
    params.setAttr("iterations", "int", 20)
    params.setAttr("mapping", "filename", "map.map")
    params.setAttr("orientation", "string", "normal")
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", 6.0)
    params.setAttr("weights", "filename", "weights.npy")

    params.setAttr("channel0_cal", "filename", "calib.npy")
    params.setAttr("channel1_cal", "filename", "calib.npy")

    params.setAttr("channel0_ext", "string", "_c1.dax")
    params.setAttr("channel1_ext", "string", "_c2.dax")

    params.setAttr("channel0_offset", "int", 0)
    params.setAttr("channel1_offset", "int", 0)

    params.setAttr("spline0", "filename", "c1_psf.spline")
    params.setAttr("spline1", "filename", "c2_psf.spline")

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.0")
    params.setAttr("max_z", "float", str(settings.spline_z_range*0.001 + 1.0))
    params.setAttr("min_z", "float", str(-settings.spline_z_range*0.001 - 1.0))

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    return params
    

# Create parameters file for analysis.
#
print("Creating XML file.")
params = testingParameters()
params.toXMLFile("multiplane.xml")

# Create localization on a grid file.
#
print("Creating gridded localization.")
sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
subprocess.call(["python", sim_path + "emitters_on_grid.py",
                 "--bin", "grid_list.bin",
                 "--nx", "14",
                 "--ny", "9",
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
                 "--bin", "psf_list.bin",
                 "--nx", "6",
                 "--ny", "3",
                 "--spacing", "40"])

# Create localization files for PSF measurement.
#
i3_locs = readinsight3.loadI3File("psf_list.bin")
for i, z_plane in enumerate(settings.z_planes):
    cx = settings.mappings["0_" + str(i) + "_x"]
    cy = settings.mappings["0_" + str(i) + "_y"]
    i3_temp = i3_locs.copy()
    xi = i3_temp["x"]
    yi = i3_temp["y"]
    xf = cx[0] + cx[1] * xi + cx[2] * yi
    yf = cy[0] + cy[1] * xi + cy[2] * yi
    i3dtype.posSet(i3_temp, "x", xf)
    i3dtype.posSet(i3_temp, "y", yf)
    i3dtype.posSet(i3_temp, "z", z_plane)
    with writeinsight3.I3Writer("c" + str(i+1) + "_psf.bin") as i3w:
        i3w.addMolecules(i3_temp)

# Create sCMOS camera calibration files.
#
numpy.save("calib.npy", [numpy.zeros((settings.x_size, settings.y_size)) + settings.camera_offset,
                         numpy.ones((settings.x_size, settings.y_size)) * settings.camera_variance,
                         numpy.ones((settings.x_size, settings.y_size)) * settings.camera_gain])

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

# Create simulated data for PSF measurements.
#
bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = 10)
cam_f = lambda s, x, y, i3 : camera.SCMOS(s, x, y, i3, "calib.npy")
drift_f = lambda s, x, y, i3 : drift.DriftFromFile(s, x, y, i3, "drift.txt")
pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, 20000.0)
psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, [])

sim = simulate.Simulate(background_factory = bg_f,
                        camera_factory = cam_f,
                        drift_factory = drift_f,
                        photophysics_factory = pp_f,
                        psf_factory = psf_f,
                        x_size = settings.x_size,
                        y_size = settings.y_size)

for i in range(len(settings.z_planes)):
    sim.simulate("c" + str(i+1) + "_zcal.dax",
                 "c" + str(i+1) + "_psf.bin",
                 dz.size)

# Create mapping file.
with open("map.map", 'wb') as fp:
    pickle.dump(settings.mappings, fp)
    
# Measure the PSF.
#
print("Measuring PSFs.")
multiplane_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/multi_plane/"
spliner_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/spliner/"
for i in range(len(settings.z_planes)):
    subprocess.call(["python", multiplane_path + "psf_zstack.py",
                     "--movie", "c" + str(i+1) + "_zcal.dax",
                     "--bin", "c" + str(i+1) + "_psf.bin",
                     "--zstack", "c" + str(i+1) + "_zstack",
                     "--scmos_cal", "calib.npy",
                     "--aoi_size", "12"])

# PSFs are independently normalized.
#
if settings.independent_heights:
    for i in range(len(settings.z_planes)):
        subprocess.call(["python", multiplane_path + "measure_psf.py",
                         "--zstack", "c" + str(i+1) + "_zstack.npy",
                         "--zoffsets", "z_offset.txt",
                         "--psf_name", "c" + str(i+1) + "_psf_normed.psf",
                         "--z_range", str(settings.spline_z_range),
                         "--normalize", "True"])

# PSFs are normalized to each other.
#
else:
    for i in range(len(settings.z_planes)):
        subprocess.call(["python", multiplane_path + "measure_psf.py",
                         "--zstack", "c" + str(i+1) + "_zstack.npy",
                         "--zoffsets", "z_offset.txt",
                         "--psf_name", "c" + str(i+1) + "_psf.psf",
                         "--z_range", str(settings.spline_z_range)])

    norm_args = ["python", multiplane_path + "normalize_psfs.py",
                 "--psfs", "c1_psf.psf"]
    for i in range(len(settings.z_planes)-1):
        norm_args.append("c" + str(i+2) + "_psf.psf")
    subprocess.call(norm_args)

# Measure the Spline.
#
print("Measuring Spline.")
#for i in range(len(settings.z_planes)):
#    subprocess.call(["python", spliner_path + "psf_to_spline.py",
#                     "--psf", "c" + str(i+1) + "_psf_normed.psf",
#                     "--spline", "c" + str(i+1) + "_psf.spline",
#                     "--spline_size", str(settings.spline_size)])

# Calculate Cramer-Rao weighting.
#
print("Calculating weights.")
subprocess.call(["python", multiplane_path + "plane_weighting.py",
                 "--background", str(settings.photons[0][0]),
                 "--photons", str(settings.photons[0][1]),
                 "--output", "weights.npy",
                 "--xml", "multiplane.xml",
                 "--no_plots", "True"])

