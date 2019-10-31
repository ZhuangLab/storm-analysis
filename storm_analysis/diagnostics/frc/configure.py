#!/usr/bin/env python
"""
Configure folder for FRC testing.

Hazen 01/18
"""
import storm_analysis.simulator.emitters_in_clusters as emittersInClusters

import storm_analysis.diagnostics.frc.settings as settings


def configure():
    # Create localizations in clusters file.
    #
    print("Creating clustered localizations file.")
    emittersInClusters.emittersInClusters("clusters_list.hdf5",
                                          50,
                                          2000,
                                          1.0,
                                          sx = settings.x_size,
                                          sy = settings.y_size)


if (__name__ == "__main__"):
    configure()
