#!/usr/bin/env python
"""
Configure folder for Multiplane PSF measurement.

Hazen 09/18
"""
import argparse
import inspect
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
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
import storm_analysis.simulator.emitters_on_grid as emittersOnGrid
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate


camera_gain = 1.0
camera_offset = 100.0
camera_variance = 1.0
independent_heights = 0
iterations = 20
pixel_size = 100.0
spline_z_range = 0.75
x_size = 300
y_size = 200
z_planes = [-0.250, 0.250]


# A mapping with a small offset.
mappings = {"0_0_x" : numpy.array([0.0, 1.0, 0.0]),
            "0_0_y" : numpy.array([0.0, 0.0, 1.0]),
            "0_1_x" : numpy.array([2.0, 1.0, 0.0]),
            "0_1_y" : numpy.array([5.0, 0.0, 1.0]),
            "1_0_x" : numpy.array([-2.0, 1.0, 0.0]),
            "1_0_y" : numpy.array([-5.0, 0.0, 1.0])}


def makeCMOSCalibration():
    # Create a CMOS calibration file.
    numpy.save("calib.npy", [numpy.zeros((y_size, x_size)) + camera_offset,
                             numpy.ones((y_size, x_size)) * camera_variance,
                             numpy.ones((y_size, x_size)) * camera_gain,
                             numpy.ones((y_size, x_size)),
                             2])


def makeMapping():
    # Create mapping file.
    with open("map.map", 'wb') as fp:
        pickle.dump(mappings, fp)

        
def makeSampleData():
    # Create sample bead data for PSF measurement.
    #

    # Create sparser grid for PSF measurement.
    #
    print("Creating data for PSF measurement.")
    emittersOnGrid.emittersOnGrid("psf_locs.hdf5", 6, 3, 1.5, 40, 0.0, 0.0)

    # Create localization files for PSF measurement.
    #
    locs = saH5Py.loadLocalizations("psf_locs.hdf5")

    for i, z_offset in enumerate(z_planes):
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
        locs_temp["z"][:] = z_offset
        
        saH5Py.saveLocalizations("c" + str(i+1) + "_psf.hdf5", locs_temp)

    # Create drift file, this is used to displace the localizations in the
    # PSF measurement movie.
    #
    dz = numpy.arange(-spline_z_range, spline_z_range + 0.001, 0.01)
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
    bg_f = lambda s, x, y, h5 : background.UniformBackground(s, x, y, h5, photons = 10)
    cam_f = lambda s, x, y, h5 : camera.SCMOS(s, x, y, h5, "calib.npy")
    drift_f = lambda s, x, y, h5 : drift.DriftFromFile(s, x, y, h5, "drift.txt")
    pp_f = lambda s, x, y, h5 : photophysics.AlwaysOn(s, x, y, h5, 20000.0)
    psf_f = lambda s, x, y, h5 : psf.PupilFunction(s, x, y, h5, pixel_size, [])

    sim = simulate.Simulate(background_factory = bg_f,
                            camera_factory = cam_f,
                            drift_factory = drift_f,
                            photophysics_factory = pp_f,
                            psf_factory = psf_f,
                            x_size = x_size,
                            y_size = y_size)

    for i in range(len(z_planes)):
        sim.simulate("c" + str(i+1) + "_zcal.dax",
                     "c" + str(i+1) + "_psf.hdf5",
                     dz.size)


def overlayImage(movie_name, locs_name, frame_number, sx = 8, sy = 8):
    """
    Create an image of a frame with the localizations overlaid.
    """
    frame = datareader.inferReader(movie_name).loadAFrame(frame_number).astype(numpy.float64)
    with saH5Py.SAH5Py(locs_name) as h5:
        locs = h5.getLocalizationsInFrame(0)

    frame = frame - numpy.min(frame)
    frame = frame/numpy.max(frame)
    
    fig = pyplot.figure(figsize = (sx, sy))
    ax = fig.add_subplot(1,1,1)
    ax.imshow(frame, interpolation = 'nearest', cmap = "gray")
    for i in range(locs["x"].size):
        width = 10
        height = 10
        if "xsigma" in locs:
            width = height = 5.0*locs["xsigma"][i]
        if "ysigma" in locs:
            height = 5.0*locs["ysigma"][i]
        ellipse = patches.Ellipse((locs["x"][i], locs["y"][i]), width, height, facecolor='none', edgecolor='g', linewidth = 2)
        ax.add_artist(ellipse)
        
    #ax.scatter(locs["x"], locs["y"], s = 200,
    ax.set_title("Overlay Image")

    pyplot.show()
    

def sCMOSSingleFrameXML():
    """
    Create a single plane sCMOS analysis file.
    """
    params = parameters.ParametersSCMOS()

    params.setAttr("max_frame", "int", 101)    
    params.setAttr("start_frame", "int", 100)    

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

    params.toXMLFile("scmos_single_plane.xml")


if (__name__ == "__main__"):
    makeCMOSCalibration()
    makeMapping()
    makeSampleData()
    sCMOSSinglePlaneXML()
    
