#!/usr/bin/env python
"""
Configure folder for Multiplane PSF mapping.

Hazen 09/18
"""
import inspect
import numpy
import os
import pickle

import storm_analysis
import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.emitters_uniform_random as emittersUniformRandom
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate


camera_gain = 1.0
camera_offset = 100.0
camera_variance = 1.0
density = 0.0005
independent_heights = 0
iterations = 20
margin = 5
pixel_size = 100.0
spline_z_range = 0.75
x_size = 300
y_size = 200
z_planes = [-0.250, 0.250]


def makeCMOSCalibration():
    # Create a CMOS calibration file.
    numpy.save("calib.npy", [numpy.zeros((y_size, x_size)) + camera_offset,
                             numpy.ones((y_size, x_size)) * camera_variance,
                             numpy.ones((y_size, x_size)) * camera_gain,
                             numpy.ones((y_size, x_size)),
                             2])

        
def makeSampleData(mappings = None):
    # Create sample bead data for mapping measurement.
    #

    # Create randomly located localizations file (for STORM movies).
    #
    print("Creating random localizations.")
    emittersUniformRandom.emittersUniformRandom("random.hdf5", density, margin, x_size, y_size, 0.0)

    # Create mapping, if not specified.
    #
    if mappings is None:
        mappings = {"0_0_x" : numpy.array([0.0, 1.0, 0.0]),
                    "0_0_y" : numpy.array([0.0, 0.0, 1.0]),
                    "0_1_x" : numpy.array([2.0, 1.0, 0.0]),
                    "0_1_y" : numpy.array([5.0, 0.0, 1.0]),
                    "1_0_x" : numpy.array([-2.0, 1.0, 0.0]),
                    "1_0_y" : numpy.array([-5.0, 0.0, 1.0])}

    # Figure out number of planes in the mapping.
    #
    n_planes = 0
    for elt in mappings:
        [i, j] = map(int, elt.split("_")[:2])
        if (i > n_planes):
            n_planes = i

    n_planes += 1
    print(n_planes)
        
    # Create localization files for PSF measurement.
    #
    locs = saH5Py.loadLocalizations("random.hdf5")

    for i in range(n_planes):
        cx = mappings["0_" + str(i) + "_x"]
        cy = mappings["0_" + str(i) + "_y"]
        locs_temp = {"x" : locs["x"].copy(),
                     "y" : locs["y"].copy(),
                     "z" : locs["z"].copy()}
        xi = locs_temp["x"]
        yi = locs_temp["y"]
        xf = cx[0] + cx[1] * xi + cx[2] * yi
        yf = cy[0] + cy[1] * xi + cy[2] * yi
        locs_temp["x"] = xf
        locs_temp["y"] = yf
        
        saH5Py.saveLocalizations("c" + str(i+1) + "_map.hdf5", locs_temp)

    # Create simulated data for PSF measurements.
    #
    bg_f = lambda s, x, y, h5 : background.UniformBackground(s, x, y, h5, photons = 10)
    cam_f = lambda s, x, y, h5 : camera.SCMOS(s, x, y, h5, "calib.npy")
    pp_f = lambda s, x, y, h5 : photophysics.AlwaysOn(s, x, y, h5, 10000.0)
    psf_f = lambda s, x, y, h5 : psf.PupilFunction(s, x, y, h5, pixel_size, [])

    sim = simulate.Simulate(background_factory = bg_f,
                            camera_factory = cam_f,
                            photophysics_factory = pp_f,
                            psf_factory = psf_f,
                            x_size = x_size,
                            y_size = y_size)

    for i in range(n_planes):
        sim.simulate("c" + str(i+1) + "_map.dax", "c" + str(i+1) + "_map.hdf5", 2)


def sCMOSXML():
    """
    Create an sCMOS analysis file.
    """
    params = parameters.ParametersSCMOS()

    params.setAttr("max_frame", "int", 1)
    params.setAttr("start_frame", "int", 0)

    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("camera_calibration", "filename", "calib.npy")
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("foreground_sigma", "float", 1.5)
    params.setAttr("iterations", "int", 1)
    params.setAttr("model", "string", "2dfixed")
    params.setAttr("pixel_size", "float", pixel_size)
    params.setAttr("roi_size", "int", 9)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", 20.0)

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.0")

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    # Z fitting.
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

    params.toXMLFile("scmos.xml")


if (__name__ == "__main__"):
    makeCMOSCalibration()
    makeSampleData()
    sCMOSXML()
    
