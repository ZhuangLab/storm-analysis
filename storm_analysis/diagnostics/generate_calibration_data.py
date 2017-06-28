#!/usr/bin/env python
"""
Generate calibration data.

This requires a list of emitter locations, like that 
created by simulator.emitters_on_grid.

/path/to/simulator/emitters_on_grid.py --bin emitters.bin --nx 6 --ny 4 --spacing 20

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


x_size = 300
y_size = 310
#z_planes = [0.0]
#z_planes = [-250.0, -250.0]
z_planes = [-250.0, 250]
#z_planes = [-750.0, -250.0, 250, 750.0]
z_range = 750.0

# Load emitter locations.
i3_locs = readinsight3.loadI3File("emitters.bin")

#
# Create channel to channel mapping file.
#
# Add offset in x,y to make this slightly more realistic.
#
mappings = {}
mappings["0_0_x"] = numpy.array([0.0, 1.0, 0.0])
mappings["0_0_y"] = numpy.array([0.0, 0.0, 1.0])
for i in range(1, len(z_planes)):
    mappings["0_" + str(i) + "_x"] = numpy.array([5.0*i, 1.0, 0.0])
    mappings["0_" + str(i) + "_y"] = numpy.array([8.0*i, 0.0, 1.0])
    #mappings["0_" + str(i) + "_y"] = numpy.array([y_size, 0.0, -1.0])
    
    mappings[str(i) + "_0_x"] = numpy.array([-5.0*i, 1.0, 0.0])
    mappings[str(i) + "_0_y"] = numpy.array([-8.0*i, 0.0, 1.0])
    #mappings[str(i) + "_0_y"] = numpy.array([y_size, 0.0, -1.0])

with open("map.map", 'wb') as fp:
    pickle.dump(mappings, fp)


# Create bin files for each plane.
for i, z_plane in enumerate(z_planes):
    cx = mappings["0_" + str(i) + "_x"]
    cy = mappings["0_" + str(i) + "_y"]
    i3_temp = i3_locs.copy()
    xi = i3_temp["x"]
    yi = i3_temp["y"]
    xf = cx[0] + cx[1] * xi + cx[2] * yi
    yf = cy[0] + cy[1] * xi + cy[2] * yi
    i3dtype.posSet(i3_temp, "x", xf)
    i3dtype.posSet(i3_temp, "y", yf)
    i3dtype.posSet(i3_temp, "z", z_plane)
    with writeinsight3.I3Writer("cam_" + str(i) + ".bin") as i3w:
        i3w.addMolecules(i3_temp)


# Create drift file, this is used to displace the
# localizations in the calibration movie.
dz = numpy.arange(-z_range, z_range + 5.0, 10.0)
drift_data = numpy.zeros((dz.size, 3))
drift_data[:,2] = dz
numpy.savetxt("drift.txt", drift_data)

# Also create the z-offset file.
z_offset = numpy.zeros((dz.size, 2))
z_offset[:,1] = dz
numpy.savetxt("z_offset.txt", z_offset)

# Create sCMOS calibration files.
offset = numpy.zeros((x_size, y_size)) + 10.0
gain = numpy.ones(offset.shape)
variance = numpy.ones(offset.shape)
for i in range(len(z_planes)):
    numpy.save("cam_cal_c" + str(i), [offset, variance, gain])

# Create simulator object.
bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = 10)
cam_f = lambda s, x, y, i3 : camera.SCMOS(s, x, y, i3, 0.0, "cam_cal_c0.npy")
drift_f = lambda s, x, y, i3 : drift.DriftFromFile(s, x, y, i3, "drift.txt")
pp_f = lambda s, x, y, i3 : photophysics.AlwaysOn(s, x, y, i3, 20000.0)

if(len(z_planes)>1):
    psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, 100.0, [])
else:
    psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, 100.0, [[1.3, 2, 2]])

sim = simulate.Simulate(background_factory = bg_f,
                        camera_factory = cam_f,
                        drift_factory = drift_f,
                        photophysics_factory = pp_f,
                        psf_factory = psf_f,
                        x_size = x_size,
                        y_size = y_size)
                        
for i in range(len(z_planes)):
    sim.simulate("zcal_c" + str(i) + ".dax",
                 "cam_" + str(i) + ".bin",
                 dz.size)

