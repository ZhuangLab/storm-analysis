#!/usr/bin/env python
"""
That which is common to all of the spliner XYZ.configure diagnostic modules.

Hazen 11/19
"""
import numpy
import pickle

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
import storm_analysis.spliner.measure_psf_utils as measurePSFUtils
import storm_analysis.spliner.psf_to_spline as psfToSpline

import storm_analysis.diagnostics.spliner.settings as settings


def configure(settings, use_dh, no_splines):
    
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

    if use_dh:
        psf_f = lambda s, x, y, i3 : psf.DHPSF(s, x, y, i3, 100.0, z_range = settings.spline_z_range)
    else:
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
    psf_name = "psf.psf"
    measurePSFBeads.measurePSFBeads("spline.dax",
                                    "z_offset.txt",
                                    "beads.txt",
                                    psf_name,
                                    aoi_size = int(settings.spline_size + 1),
                                    pixel_size = settings.pixel_size * 1.0e-3,
                                    z_range = settings.spline_z_range)

    # Smooth the PSF if requested.
    #
    if settings.smooth_psf:
        with open(psf_name, "rb") as fp:
            psf_data = pickle.load(fp)
        sm_psf = measurePSFUtils.smoothPSF(psf_data["psf"],
                                           xy_sigma = settings.smooth_psf_sigma,
                                           z_sigma = settings.smooth_psf_sigma)
        psf_data["psf"] = sm_psf
        psf_data["smoothed"] = True

        psf_name = "psf_smooth.psf"
        with open(psf_name, "wb") as fp:
            pickle.dump(psf_data, fp)
            
    # Measure the Spline.
    #
    if True:
        print("Measuring Spline.")
        psfToSpline.psfToSpline(psf_name, "psf.spline", settings.spline_size)
