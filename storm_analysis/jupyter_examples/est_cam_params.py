#!/usr/bin/env python
"""
Configure folder for camera parameter estimation Jupyter example.

Hazen 09/18
"""
import inspect
import os

import storm_analysis
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.emitters_uniform_random as emittersUniformRandom
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

bg = 20
margin = 1
n_frames = 200
pixel_size = 100.0
signal = 2000 # in photons.
x_size = 300
y_size = 250

def createLocalizations():
    emittersUniformRandom.emittersUniformRandom("sim_locs.hdf5", 0.1, margin, x_size, y_size, 0.0) 

def createMovie(gain = 1.0, offset = 100.0):
    bg_f = lambda s, x, y, h5 : background.GaussianBackground(s, x, y, h5, photons = bg)
    cam_f = lambda s, x, y, h5 : camera.Ideal(s, x, y, h5, offset, gain = gain)
    pp_f = lambda s, x, y, h5 : photophysics.SimpleSTORM(s, x, y, h5, photons = signal)
    psf_f = lambda s, x, y, h5 : psf.GaussianPSF(s, x, y, h5, pixel_size)

    sim = simulate.Simulate(background_factory = bg_f,
                            camera_factory = cam_f,
                            photophysics_factory = pp_f,
                            psf_factory = psf_f,
                            x_size = x_size,
                            y_size = y_size)
    
    sim.simulate("test.tif", "sim_locs.hdf5", n_frames)

def configure():
    print("Creating ground truth localizations.")
    createLocalizations()
    print("Creating movie.")
    createMovie()


if (__name__ == "__main__"):
    configure()
