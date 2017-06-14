#!/usr/bin/env python
"""
Generate calibration data.

This requires a list of emitter locations, like that 
created by simulator.emitters_on_grid.

/path/to/simulator/emitters_on_grid.py --bin emitters.bin --nx 6 --ny 4 --spacing 20

It also assume that:
/path/to/simulator/generate_calibration_data.py

has been run in this directory.

Run in the directory with the emitters.bin file.
"""

import numpy
import pickle

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate


frames = 100
x_size = 300
y_size = 200
z_planes = [0.0]
#z_planes = [-250.0, 250]
#z_planes = [-750.0, -250.0, 250, 750.0]
z_value = 0.0

# Load emitter locations.
i3_locs = readinsight3.loadI3File("emitters.bin")

# Create bin files for each plane.
for i, z_plane in enumerate(z_planes):
    i3dtype.posSet(i3_locs, "z", z_plane + z_value)
    with writeinsight3.I3Writer("sim_input_c" + str(i) + ".bin") as i3w:
        i3w.addMolecules(i3_locs)

# Create simulator object.
bg_photons = int(100.0/float(len(z_planes)))
signal = 6000.0/float(len(z_planes))
    
bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg_photons)
cam_f = lambda s, x, y, i3 : camera.SCMOS(s, x, y, i3, 0.0, "cam_cal_c0.npy")
pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, signal)

if(len(z_planes)>1):
    psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, 100.0, [])
else:
    psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, 100.0, [[1.3, 2, 2]])

sim = simulate.Simulate(background_factory = bg_f,
                        camera_factory = cam_f,
                        photophysics_factory = pp_f,
                        psf_factory = psf_f,
                        x_size = x_size,
                        y_size = y_size)
                        
for i in range(len(z_planes)):
    sim.simulate("test_c" + str(i) + ".dax",
                 "sim_input_c" + str(i) + ".bin",
                 frames)

