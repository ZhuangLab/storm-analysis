#!/usr/bin/env python
"""
Configure folder for testing Cramer-Rao bounds calculations. Ideally 
this should be a unit test, but calculating splines takes too long and
they are a bit big to include in a project, so such is life.

Hazen 10/17
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

import settings


#parser = argparse.ArgumentParser(description = 'Cramer-Rao diagnostics configuration.')
#args = parser.parse_args()

def configure():

    # Create PF for pupil function.
    #
    pupilfn_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/pupilfn/"
    print("Creating pupil function.")
    subprocess.call(["python", pupilfn_path + "make_pupil_fn.py",
                     "--filename", "pupilfn.pfn",
                     "--size", str(settings.spline_size),
                     "--pixel-size", str(settings.pixel_size),
                     "--zmn", str(settings.zmn),
                     "--z-offset", str(settings.z_offset)])
    
    # Create PSF using pupil functions directly.
    #
    if False:
        psf_fft_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/psf_fft/"
        print("Creating (theoritical) psf.")
        subprocess.call(["python", psf_fft_path + "make_psf_from_pf.py",
                         "--filename", "psf_fft.psf",
                         "--size", str(settings.spline_size),
                         "--pixel-size", str(settings.pixel_size),
                         "--zrange", str(settings.psf_fft_z_range),
                         "--zstep", str(settings.psf_fft_z_step)])

        exit()

    # Localizations on a sparse parse grid for PSF
    # measurement for Spliner and PSF FFT.
    #
    print("Creating data for PSF measurement.")
    sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
    subprocess.call(["python", sim_path + "emitters_on_grid.py",
                     "--bin", "sparse_list.hdf5",
                     "--nx", "6",
                     "--ny", "3",
                     "--spacing", "40",
                     "--zoffset", str(settings.z_offset)])

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
                        
    sim.simulate("psf.dax", "sparse_list.hdf5", dz.size)

    # Create spline for Spliner
    #

    # Measure the PSF for Spliner
    #
    print("Measuring PSF.")
    spliner_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/spliner/"
    subprocess.call(["python", spliner_path + "measure_psf_beads.py",
                     "--movie", "psf.dax",
                     "--zoffset", "z_offset.txt",
                     "--aoi_size", str(int(settings.spline_size/2)+1),
                     "--beads", "beads.txt",
                     "--psf", "psf_spliner.psf"])

    # Measure the Spline.
    #
    
    # This is slow, sometimes you don't want to do it.
    if True:
        print("Measuring Spline.")
        subprocess.call(["python", spliner_path + "psf_to_spline.py",
                         "--psf", "psf_spliner.psf",
                         "--spline", "psf.spline",
                         "--spline_size", str(settings.spline_size)])
        

    # Create measured PSF for PSF FFT.
    #

    # Measure the PSF using spliner/measure_psf_beads.py
    #
    print("Measuring PSF.")
    subprocess.call(["python", spliner_path + "measure_psf_beads.py",
                     "--movie", "psf.dax",
                     "--zoffset", "z_offset.txt",
                     "--aoi_size", str(int(settings.spline_size/2)+1),
                     "--beads", "beads.txt",
                     "--psf", "psf_fft_2x.psf",
                     "--zrange", str(settings.psf_fft_z_range),
                     "--zstep", str(settings.psf_fft_z_step)])

    # Downsample by 2x for use by PSF FFT.
    #
    psf_fft_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/psf_fft/"
    print("Creating downsampled psf.")
    subprocess.call(["python", psf_fft_path + "downsample_psf.py",
                     "--spliner_psf", "psf_fft_2x.psf",
                     "--psf", "psf_fft.psf",
                     "--pixel-size", str(settings.pixel_size)])


if (__name__ == "__main__"):
    configure()
    
