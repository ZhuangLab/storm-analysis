#!/usr/bin/env python
"""
Configure folder for drift correction testing.

Note: This test takes approximately 20-30 minutes.

Hazen 01/18
"""
import numpy

import storm_analysis.sa_library.parameters as parameters

import storm_analysis.simulator.emitters_in_clusters as emittersInClusters
import storm_analysis.simulator.emitters_on_lines as emittersOnLines

import storm_analysis.diagnostics.drift_correction.settings as settings


def testingParameters():
    """
    Create a 3D-DAOSTORM parameters object.
    """
    params = parameters.ParametersDAO()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)
    
    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("camera_gain", "float", settings.camera_gain)
    params.setAttr("camera_offset", "float", settings.camera_offset)
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("foreground_sigma", "float", 1.0)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("model", "string", settings.model)
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("roi_size", "int", 9)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", 6.0)

    # Do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.5")

    # Do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 1)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 1)

    return params


def configure():
    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters()
    params.toXMLFile("dao.xml")

    # Create localizations in clusters file.
    #
    print("Creating clustered localizations file.")
    emittersInClusters.emittersInClusters("clusters_list.hdf5",
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
    
