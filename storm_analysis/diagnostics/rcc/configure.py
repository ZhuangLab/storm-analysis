#!/usr/bin/env python
"""
Configure folder for RCC drift correction testing.

Note: This test takes approximately 20-30 minutes.

Hazen 09/18
"""
import numpy

import storm_analysis.simulator.emitters_in_clusters as emittersInClusters
import storm_analysis.simulator.emitters_on_lines as emittersOnLines

import storm_analysis.diagnostics.rcc.settings as settings

    
def configure():
    
    # Create localizations in clusters file.
    #
    print("Creating clustered localizations file.")
    emittersInClusters.emittersInClusters("cluster_list.hdf5",
                                          50,
                                          200,
                                          1.0,
                                          sx = settings.x_size,
                                          sy = settings.y_size)

    # Create localizations on lines file.
    #
    print("Creating lines localizations file.")
    emittersOnLines.emittersOnLines("lines_list.hdf5",
                                    50,
                                    100000,
                                    sx = settings.x_size,
                                    sy = settings.y_size)

    # Create drift file. This is used in the simulations to displace
    # the localizations.
    #
    dx = settings.x_drift/settings.n_frames
    drift_x = numpy.arange(0.0, settings.x_drift + 0.5 * dx, dx)
    
    dy = settings.y_drift/settings.n_frames
    drift_y = numpy.arange(0.0, settings.y_drift + 0.5 * dy, dy)

    dz = settings.z_drift/settings.n_frames
    drift_z = numpy.arange(0.0, settings.z_drift + 0.5 * dz, dz)

    drift_data = numpy.zeros((drift_x.size, 3))
    drift_data[:,0] = drift_x
    drift_data[:,1] = drift_y
    
    numpy.savetxt("drift_xy.txt", drift_data)
    
    drift_data[:,2] = drift_z
    numpy.savetxt("drift_xyz.txt", drift_data)


if (__name__ == "__main__"):
    configure()
    
