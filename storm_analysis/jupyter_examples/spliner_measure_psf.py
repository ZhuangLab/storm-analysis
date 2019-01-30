#!/usr/bin/env python
"""
Configure folder for Spliner PSF measurement.

Hazen 03/18
"""
import inspect
import numpy
import os

import storm_analysis
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

bg = 20
camera_gain = 1.0
camera_offset = 100.0
iterations = 20
margin = 1
model = "3d"
nx = 12
ny = 10
pixel_size = 100.0
signal = 4000 # in photons.
tolerance = 0.3
x_size = 256
y_size = 256
z_range = 0.6


def createLocalizations():
    """
    Create a 'realistic' localizations file with small relative z offsets.

    We'll also use these to create the bead locations text file.
    """
    x = numpy.array([30.3, 43.1, 133.2, 111.7, 200.2])
    y = numpy.array([205.4, 56.2, 176.4, 93.6, 134.0])
    z = numpy.array([-0.05, -0.02, 0.0, 0.02, 0.05])
    sx = numpy.ones(x.size)
    sy = numpy.ones(x.size)
    
    with saH5Py.SAH5Py("spliner_measure_psf.hdf5", is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(x_size, y_size, 1, "")
        h5.setPixelSize(pixel_size)
        h5.addLocalizations({"x" : x, "y" : y, "z" : z, "xsigma" : sx, "ysigma" : sy}, 0)

    blocs = numpy.zeros((x.size, 2))
    blocs[:,0] = x + numpy.random.uniform(low = -0.2, high = 0.2, size = x.size)
    blocs[:,1] = y + numpy.random.uniform(low = -0.2, high = 0.2, size = x.size)
    numpy.savetxt("bead_locs.txt", blocs)
    
    
def createMovie(movie_name, use_dh = False):
    
    # Create drift file, this is used to displace the localizations in the
    # PSF measurement movie.
    #
    dz = numpy.arange(-z_range, z_range + 0.001, 0.01)
    drift_data = numpy.zeros((dz.size, 3))
    drift_data[:,2] = dz
    numpy.savetxt("drift.txt", drift_data)

    # Also create the z-offset file.
    #
    z_offset = numpy.ones((dz.size, 2))
    z_offset[:,1] = dz
    numpy.savetxt("z_offsets.txt", z_offset)
    
    bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
    cam_f = lambda s, x, y, i3 : camera.Ideal(s, x, y, i3, camera_offset)
    drift_f = lambda s, x, y, i3 : drift.DriftFromFile(s, x, y, i3, "drift.txt")
    pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, signal)
    if use_dh:
        psf_f = lambda s, x, y, i3 : psf.DHPSF(s, x, y, i3, pixel_size)
    else:
        psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, pixel_size, [[1.3, 2, 2]])

    sim = simulate.Simulate(background_factory = bg_f,
                            camera_factory = cam_f,
                            drift_factory = drift_f,
                            photophysics_factory = pp_f,
                            psf_factory = psf_f,
                            x_size = x_size,
                            y_size = y_size)
    
    sim.simulate(movie_name, "spliner_measure_psf.hdf5", dz.size)

    
def configure():
    print("Creating bead localizations.")
    createLocalizations()
    print("Creating measurement movie.")
    createMovie("spliner_measure.tif")


if (__name__ == "__main__"):
    configure()
