#!/usr/bin/env python
"""
Configure folder for SLURM testing.

Hazen 08/18
"""
import numpy
import os
import pickle
import shutil

import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.emitters_on_grid as emittersOnGrid
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.multi_plane.plane_weighting as planeWeighting

import storm_analysis.pupilfn.make_pupil_fn as makePupilFn

import storm_analysis.diagnostics.slurm.settings as settings


def testingParameters():
    """
    Create a Multiplane parameters object.
    """
    params = parameters.ParametersMultiplaneArb()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)
    
    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("find_max_radius", "int", 2)
    params.setAttr("independent_heights", "int", settings.independent_heights)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("mapping", "filename", "map.map")
    params.setAttr("no_fitting", "int", 0)
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", 6.0)
    params.setAttr("weights", "filename", "weights.npy")
    params.setAttr("z_value", "float-array", settings.z_value)

    params.setAttr("channel0_cal", "filename", "calib.npy")
    params.setAttr("channel1_cal", "filename", "calib.npy")

    params.setAttr("channel0_ext", "string", "_c1.dax")
    params.setAttr("channel1_ext", "string", "_c2.dax")

    params.setAttr("channel0_offset", "int", 0)
    params.setAttr("channel1_offset", "int", 0)

    params.setAttr("pupilfn0", "filename", "c1_pupilfn.pfn")
    params.setAttr("pupilfn1", "filename", "c2_pupilfn.pfn")
        
    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.0")
    
    params.setAttr("max_z", "float", str(settings.pupilfn_z_range))
    params.setAttr("min_z", "float", str(-settings.pupilfn_z_range))

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    return params


def configure():
    
    # Create directory, if necessary.
    if not os.path.exists(settings.wdir):
        os.makedirs(settings.wdir)
    
    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters()
    params.toXMLFile("multiplane.xml")

    # Create localization on a grid file.
    #
    print("Creating gridded localization.")
    emittersOnGrid.emittersOnGrid("grid_list.hdf5",
                                  settings.nx,
                                  settings.ny,
                                  1.5,
                                  20,
                                  settings.test_z_range,
                                  settings.test_z_offset)

    # Create sCMOS camera calibration files.
    #
    numpy.save(os.path.join(settings.wdir, "calib.npy"),
               [numpy.zeros((settings.y_size, settings.x_size)) + settings.camera_offset,
                numpy.ones((settings.y_size, settings.x_size)) * settings.camera_variance,
                numpy.ones((settings.y_size, settings.x_size)) * settings.camera_gain,
                numpy.ones((settings.y_size, settings.x_size)),
                2])
    shutil.copyfile(os.path.join(settings.wdir, "calib.npy"), "calib.npy")

    # Create mapping file.
    with open(os.path.join(settings.wdir, "map.map"), 'wb') as fp:
        pickle.dump(settings.mappings, fp)
    shutil.copyfile(os.path.join(settings.wdir, "map.map"), "map.map")

    # Create pupil functions for 'pupilfn'.
    print("Creating pupil functions.")
    for i in range(len(settings.z_planes)):
        fname = "c" + str(i+1) + "_pupilfn.pfn"
        makePupilFn.makePupilFunction(os.path.join(settings.wdir, fname),
                                      settings.psf_size,
                                      settings.pixel_size * 1.0e-3,
                                      settings.pupil_fn,
                                      z_offset = -settings.z_planes[i])

        shutil.copyfile(os.path.join(settings.wdir, fname), fname)

    # Calculate Cramer-Rao weighting.
    #
    print("Calculating weights.")
    planeWeighting.runPlaneWeighting("multiplane.xml",
                                     os.path.join(settings.wdir, "weights.npy"),
                                     [settings.photons[0]],
                                     settings.photons[1],
                                     no_plots = True)


if (__name__ == "__main__"):
    configure()
