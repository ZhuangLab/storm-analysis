#!/usr/bin/env python
"""
Configure folder for testing Cramer-Rao bounds calculations. Ideally 
this should be a unit test, but calculating splines takes too long and
they are a bit big to include in a project, so such is life.

Hazen 10/17
"""
import argparse
import numpy
import os

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

import storm_analysis.pupilfn.make_pupil_fn as makePupilFn

import storm_analysis.psf_fft.make_psf_from_pf as makePSFFromPF

import storm_analysis.spliner.measure_psf_beads as measurePSFBeads
import storm_analysis.spliner.psf_to_spline as psfToSpline

import storm_analysis.diagnostics.cramer_rao.settings as settings


#parser = argparse.ArgumentParser(description = 'Cramer-Rao diagnostics configuration.')
#args = parser.parse_args()

def configure():

    # Create PF for pupil function.
    #
    print("Creating pupil function.")
    pf_size = 2*(settings.spline_size - 1)
    makePupilFn.makePupilFunction("pupilfn.pfn",
                                  pf_size,
                                  settings.pixel_size * 1.0e-3,
                                  settings.zmn,
                                  z_offset = settings.z_offset)

    # Create PSF using pupil functions directly.
    #
    if False:
        print("Creating (theoritical) psf.")
        makePSFFromPF.makePSF("psf_fft.psf",
                              settings.spline_size,
                              settings.pixel_size * 1.0e-3,
                              settings.zmn,
                              settings.psf_fft_z_range,
                              settings.psf_fft_z_step)

        exit()

    # Localizations on a sparse parse grid for PSF
    # measurement for Spliner and PSF FFT.
    #
    print("Creating data for PSF measurement.")
    emittersOnGrid.emittersOnGrid("sparse_list.hdf5",
                                  6,
                                  3,
                                  1.5,
                                  40,
                                  0.0,
                                  settings.z_offset)

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
    cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, 100.0)
    drift_f = lambda s, x, y, i3 : drift.DriftFromFile(s, x, y, i3, "drift.txt")
    pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, 20000.0)
    psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, settings.zmn, pf_size = pf_size)

    sim = simulate.Simulate(background_factory = bg_f,
                            camera_factory = cam_f,
                            drift_factory = drift_f,
                            photophysics_factory = pp_f,
                            psf_factory = psf_f,
                            x_size = settings.x_size,
                            y_size = settings.y_size)
                        
    sim.simulate("psf.dax", "sparse_list.hdf5", dz.size)

    # Create spline for Spliner
    #

    # Measure the PSF for Spliner
    #
    print("Measuring PSF.")
    psf_name = "psf_spliner.psf"
    measurePSFBeads.measurePSFBeads("psf.dax",
                                    "z_offset.txt",
                                    "beads.txt",
                                    psf_name,
                                    aoi_size = int(settings.spline_size + 1),
                                    pixel_size = settings.pixel_size * 1.0e-3)

    # Measure the Spline.
    #
    
    # This is slow, sometimes you don't want to do it.
    if True:
        print("Measuring Spline.")
        psfToSpline.psfToSpline(psf_name, "psf.spline", settings.spline_size)

    # Create measured PSF for PSF FFT.
    #

    # Measure the PSF using spliner/measure_psf_beads.py
    #
    print("Measuring PSF.")
    measurePSFBeads.measurePSFBeads("psf.dax",
                                    "z_offset.txt",
                                    "beads.txt",
                                    "psf_fft.psf",
                                    aoi_size = int(settings.spline_size - 1),
                                    pixel_size = settings.pixel_size * 1.0e-3,
                                    z_range = settings.psf_fft_z_range,
                                    z_step = settings.psf_fft_z_step)


if (__name__ == "__main__"):
    configure()
    
