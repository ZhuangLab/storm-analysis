#!/usr/bin/env python
"""
Configure folder for RCC drift correction testing.

Note: This test takes approximately 20-30 minutes.

Hazen 09/18
"""
import inspect
import numpy
import os
import shutil
import subprocess

import storm_analysis
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.diagnostics.rcc.settings as settings

    
def configure():
    
    # Create localizations in clusters file.
    #
    print("Creating clustered localizations file.")
    sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
    subprocess.call(["python", sim_path + "emitters_in_clusters.py",
                     "--bin", "clusters_list.hdf5",
                     "--dev", "1.0",
                     "--ncl", "50",
                     "--nlocs", "200",
                     "--sx", str(settings.x_size),
                     "--sy", str(settings.y_size)])

    # Create localizations on lines file.
    #
    print("Creating lines localizations file.")
    subprocess.call(["python", sim_path + "emitters_on_lines.py",
                     "--bin", "lines_list.hdf5",
                     "--nlines", "50",
                     "--nemitters", "100000",
                     "--sx", str(settings.x_size),
                     "--sy", str(settings.y_size)])

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
    
