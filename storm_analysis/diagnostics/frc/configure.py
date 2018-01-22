#!/usr/bin/env python
"""
Configure folder for FRC testing.

Hazen 01/18
"""
import inspect
import os
import subprocess

import storm_analysis

import settings


# Create localizations in clusters file.
#
print("Creating clustered localizations file.")
sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
subprocess.call(["python", sim_path + "emitters_in_clusters.py",
                 "--bin", "clusters_list.hdf5",
                 "--dev", "1.0",
                 "--ncl", "50",
                 "--nlocs", "2000",
                 "--sx", str(settings.x_size),
                 "--sy", str(settings.y_size)])
