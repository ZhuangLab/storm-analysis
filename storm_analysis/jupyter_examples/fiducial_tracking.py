#!/usr/bin/env python
"""
Configure folder for fiducials tracking.

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


camera_offset = 100.0
density = 0.0005
iterations = 1
margin = 5
n_frames = 100
pixel_size = 100.0
threshold = 10.0
x_size = 300
y_size = 200

# Constants to use to generate the drift data.
x_loc = 0.025
x_scale = 0.1
y_loc = 0.1 * x_loc
y_scale = 0.1


def makeSampleData():
    # Create sample bead data for fiducial tracking.
    #

    # Create randomly located localizations file.
    #
    print("Creating random localizations.")
    emittersUniformRandom.emittersUniformRandom("random.hdf5", density, margin, x_size, y_size, 0.0)
    
    # Create X/Y/Z drift file.
    #
    dx = numpy.random.normal(loc = x_loc, scale = x_scale, size = n_frames)
    dy = numpy.random.normal(loc = y_loc, scale = y_scale, size = n_frames)

    # Integrate dx, dy
    for i in range(1, dx.size):
        dx[i] = dx[i-1] + dx[i]
        dy[i] = dy[i-1] + dy[i]
    
    drift_data = numpy.zeros((dx.size, 3))
    drift_data[:,0] = dx
    drift_data[:,1] = dy
    numpy.savetxt("drift.txt", drift_data)    
        
    # Create simulated data for fiducial tracking.
    #
    bg_f = lambda s, x, y, h5 : background.UniformBackground(s, x, y, h5, photons = 10)
    cam_f = lambda s, x, y, h5 : camera.Ideal(s, x, y, h5, camera_offset)
    drift_f = lambda s, x, y, h5 : drift.DriftFromFile(s, x, y, h5, "drift.txt")    
    pp_f = lambda s, x, y, h5 : photophysics.SimpleSTORM(s, x, y, h5, 4000, on_time = 10.0, off_time = 1.0)
    psf_f = lambda s, x, y, h5 : psf.GaussianPSF(s, x, y, h5, pixel_size)

    sim = simulate.Simulate(background_factory = bg_f,
                            camera_factory = cam_f,
                            drift_factory = drift_f,
                            photophysics_factory = pp_f,
                            psf_factory = psf_f,
                            x_size = x_size,
                            y_size = y_size)
    
    sim.simulate("fiducials.tif", "random.hdf5", n_frames)


def daoSTORMXML():
    """
    Create a 3D-DAOSTORM parameters object.
    """
    params = parameters.ParametersDAO()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)

    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("camera_gain", "float", 1.0)
    params.setAttr("camera_offset", "float", camera_offset)
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("foreground_sigma", "float", 1.0)
    params.setAttr("iterations", "int", iterations)
    params.setAttr("model", "string", "2dfixed")
    params.setAttr("pixel_size", "float", pixel_size)
    params.setAttr("roi_size", "int", 9)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", threshold)

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0")

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

    params.toXMLFile("daostorm.xml")


if (__name__ == "__main__"):
    makeSampleData()
    daoSTORMXML()
    
