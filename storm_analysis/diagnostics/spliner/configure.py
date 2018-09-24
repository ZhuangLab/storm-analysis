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
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.diagnostics.spliner.settings as settings


def testingParameters(cal_file = None):
    """
    Create a Spliner parameters object.
    """
    params = parameters.ParametersSpliner()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)    

    params.setAttr("anscombe", "int", settings.anscombe)
    params.setAttr("background_sigma", "float", 8.0)

    if cal_file is not None:
        params.setAttr("camera_calibration", "filename", cal_file)
    else:
        params.setAttr("camera_gain", "float", settings.camera_gain)
        params.setAttr("camera_offset", "float", settings.camera_offset)

    params.setAttr("find_max_radius", "int", 5)
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
    
    # Create sparser grid for PSF measurement.
    #
    print("Creating data for PSF measurement.")
    sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
    subprocess.call(["python", sim_path + "emitters_on_grid.py",
                     "--bin", "sparse_list.hdf5",
                     "--nx", "6",
                     "--ny", "3",
                     "--spacing", "40"])
        
    if no_splines:
        return
    
    # Create beads.txt file for spline measurement.
    #
    with saH5Py.SAH5Py("sparse_list.hdf5") as h5:
        locs = h5.getLocalizations()
        numpy.savetxt("beads.txt", numpy.transpose(numpy.vstack((locs['x'], locs['y']))))

    # Create drift file, this is used to displace the localizations in the
    # PSF measurement movie.
    #
    dz = numpy.arange(-settings.spline_z_range, settings.spline_z_range + 0.001, 0.01)
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
                        
    sim.simulate("spline.dax", "sparse_list.hdf5", dz.size)

    # Measure the PSF.
    #
    print("Measuring PSF.")
    spliner_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/spliner/"
    subprocess.call(["python", spliner_path + "measure_psf_beads.py",
                     "--movie", "spline.dax",
                     "--zoffset", "z_offset.txt",
                     "--beads", "beads.txt",
                     "--psf", "psf.psf",
                     "--aoi_size", str(settings.spline_size+1)])

    # Measure the Spline.
    #
    if True:
        print("Measuring Spline.")
        subprocess.call(["python", spliner_path + "psf_to_spline.py",
                         "--psf", "psf.psf",
                         "--spline", "psf.spline",
                         "--spline_size", str(settings.spline_size)])


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description = 'Spline diagnostics configuration.')

    parser.add_argument('--no-splines', dest='no_splines', action='store_true', default = False)

    args = parser.parse_args()
    
    configure(args.no_splines)
