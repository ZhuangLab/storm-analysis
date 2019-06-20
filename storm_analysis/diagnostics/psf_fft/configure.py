#!/usr/bin/env python
"""
Configure folder for PSF FFT testing.

Hazen 10/17
"""
import numpy
import os

import storm_analysis
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.emitters_on_grid as emittersOnGrid
import storm_analysis.simulator.emitters_uniform_random as emittersUniformRandom
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.spliner.measure_psf_beads as measurePSFBeads

import storm_analysis.psf_fft.make_psf_from_pf as makePSFFromPF

import storm_analysis.diagnostics.psf_fft.settings as settings


def testingParameters(cal_file = None):
    """
    Create a PSF FFT parameters object.
    """
    params = parameters.ParametersPSFFFT()

    params.setAttr("max_frame", "int", -1) 
    params.setAttr("start_frame", "int", -1)
    
    params.setAttr("background_sigma", "float", 8.0)

    if cal_file is not None:
        params.setAttr("camera_calibration", "filename", cal_file)
    else:
        params.setAttr("camera_gain", "float", settings.camera_gain)
        params.setAttr("camera_offset", "float", settings.camera_offset)
        
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("iterations", "int", settings.iterations)
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

    # 'peak_locations' testing.
    if hasattr(settings, "peak_locations") and (settings.peak_locations is not None):
        params.setAttr("peak_locations", "filename", settings.peak_locations)

    return params


def configure(cal_file = None):
    
    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters(cal_file = cal_file)
    params.toXMLFile("psf_fft.xml")

    # Create localization on a grid file.
    #
    print("Creating gridded localization.")
    emittersOnGrid.emittersOnGrid("grid_list.hdf5",
                                  settings.nx,
                                  settings.ny,
                                  1.5,
                                  20,
                                  settings.test_z_range,
                                  settings.test_z_offset)

    # Create randomly located localizations file.
    #
    print("Creating random localization.")
    emittersUniformRandom.emittersUniformRandom("random_list.hdf5",
                                                1.0,
                                                settings.margin,
                                                settings.x_size,
                                                settings.y_size,
                                                settings.test_z_range)

    # Create sparser grid for PSF measurement.
    #
    print("Creating data for PSF measurement.")
    emittersOnGrid.emittersOnGrid("sparse_list.hdf5",
                                  6,
                                  3,
                                  1.5,
                                  40,
                                  0.0,
                                  0.0)


    if False:
    
        # Create PSF using pupil functions directly.
        #
        print("Creating (theoritical) psf.")
        makePSFFromPF.makePSF("psf.psf",
                              settings.psf_size,
                              settings.pixel_size * 1.0e-3,
                              settings.zmn,
                              settings.psf_z_range,
                              settings.z_step)

    else:

        # Create beads.txt file for PSF measurement.
        #
        with saH5Py.SAH5Py("sparse_list.hdf5") as h5:
            locs = h5.getLocalizations()
            numpy.savetxt("beads.txt", numpy.transpose(numpy.vstack((locs['x'], locs['y']))))

        # Create drift file, this is used to displace the localizations in the
        # PSF measurement movie.
        #
        dz = numpy.arange(-settings.psf_z_range, settings.psf_z_range + 0.001, 0.010)
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
    
        sim.simulate("psf.dax", "sparse_list.hdf5", dz.size)

        # Measure the PSF using spliner/measure_psf_beads.py
        #
        print("Measuring PSF.")
        measurePSFBeads.measurePSFBeads("psf.dax",
                                        "z_offset.txt",
                                        "beads.txt",
                                        "psf.psf",
                                        aoi_size = int(settings.psf_size/2)+1,
                                        pixel_size = settings.pixel_size * 1.0e-3,
                                        z_range = settings.psf_z_range,
                                        z_step = settings.z_step)


if (__name__ == "__main__"):
    configure()
    
