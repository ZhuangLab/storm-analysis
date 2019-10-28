#!/usr/bin/env python
"""
Configure folder for Multicolor testing.

Hazen 01/18
"""
import argparse
import numpy
import os
import pickle

import storm_analysis
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.emitters_on_grid as emittersOnGrid
import storm_analysis.simulator.emitters_uniform_random as emittersUniformRandom
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.micrometry.merge_maps as mergeMaps
import storm_analysis.micrometry.micrometry as micrometry

import storm_analysis.multi_plane.measure_psf as mpMeasurePSF
import storm_analysis.multi_plane.normalize_psfs as normalizePSFs
import storm_analysis.multi_plane.plane_weighting as planeWeighting
import storm_analysis.multi_plane.print_mapping as printMapping
import storm_analysis.multi_plane.psf_localizations as psfLocalizations
import storm_analysis.multi_plane.psf_zstack as psfZStack

import storm_analysis.sCMOS.scmos_analysis as scmos

import storm_analysis.spliner.psf_to_spline as psfToSpline

import storm_analysis.diagnostics.multicolor.settings as settings


def testingParametersSCMOS():
    """
    Create a sCMOS parameters object.
    """
    params = parameters.ParametersSCMOS()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)    
    
    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("camera_calibration", "filename", "calib.npy")
    params.setAttr("find_max_radius", "int", 5)
    params.setAttr("foreground_sigma", "float", 1.5)
    params.setAttr("iterations", "int", settings.iterations)
    params.setAttr("model", "string", "2dfixed")
    params.setAttr("pixel_size", "float", settings.pixel_size)
    params.setAttr("sigma", "float", 150.0/settings.pixel_size)
    params.setAttr("threshold", "float", 6.0)

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.0")

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)
    
    return params


def testingParametersMC():
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
    params.setAttr("channel2_cal", "filename", "calib.npy")
    params.setAttr("channel3_cal", "filename", "calib.npy")

    params.setAttr("channel0_ext", "string", "_c1.dax")
    params.setAttr("channel1_ext", "string", "_c2.dax")
    params.setAttr("channel2_ext", "string", "_c3.dax")
    params.setAttr("channel3_ext", "string", "_c4.dax")

    params.setAttr("channel0_offset", "int", 0)
    params.setAttr("channel1_offset", "int", 0)
    params.setAttr("channel2_offset", "int", 0)
    params.setAttr("channel3_offset", "int", 0)

    params.setAttr("spline0", "filename", "c1_psf.spline")
    params.setAttr("spline1", "filename", "c2_psf.spline")
    params.setAttr("spline2", "filename", "c3_psf.spline")
    params.setAttr("spline3", "filename", "c4_psf.spline")    

    # Do tracking (localization color analysis depends on the tracks).
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "1.0")

    params.setAttr("max_z", "float", str(0.001 * settings.psf_z_range))
    params.setAttr("min_z", "float", str(-0.001 * settings.psf_z_range))
    
    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    return params

def configure():
    # Create analysis XML files.
    #
    print("Creating XML files.")
    params = testingParametersSCMOS()
    params.toXMLFile("scmos.xml")

    params = testingParametersMC()
    params.toXMLFile("multicolor.xml")
    
    # Useful variables
    aoi_size = int(settings.psf_size/2)+1

    # Create sCMOS data and HDF5 files we'll need for the simulation.
    #
    if True:

        # Create sCMOS camera calibration files.
        #
        numpy.save("calib.npy", [numpy.zeros((settings.y_size, settings.x_size)) + settings.camera_offset,
                                 numpy.ones((settings.y_size, settings.x_size)) * settings.camera_variance,
                                 numpy.ones((settings.y_size, settings.x_size)) * settings.camera_gain,
                                 numpy.ones((settings.y_size, settings.x_size)),
                                 2])
    
        # Create localization on a grid file.
        #
        print("Creating gridded localizations.")
        emittersOnGrid.emittersOnGrid("grid_list.hdf5",
                                      settings.nx,
                                      settings.ny,
                                      1.5,
                                      20,
                                      settings.test_z_range,
                                      settings.test_z_offset)

        # Create randomly located localizations file (for STORM movies).
        #
        print("Creating random localizations.")
        emittersUniformRandom.emittersUniformRandom("random_storm.hdf5",
                                                    1.0,
                                                    settings.margin,
                                                    settings.x_size,
                                                    settings.y_size,
                                                    settings.test_z_range)

        # Create randomly located localizations file (for mapping measurement).
        #
        print("Creating random localizations.")
        emittersUniformRandom.emittersUniformRandom("random_map.hdf5",
                                                    0.0003,
                                                    settings.margin,
                                                    settings.x_size,
                                                    settings.y_size,
                                                    settings.test_z_range)

        # Create sparser grid for PSF measurement.
        #
        print("Creating data for PSF measurement.")
        emittersOnGrid.emittersOnGrid("psf_list.hdf5",
                                      6,
                                      3,
                                      1.5,
                                      40,
                                      0.0,
                                      0.0)


    ## This part makes / tests measuring the mapping.
    ##
    if True:
        print("Measuring mapping.")
    
        # Make localization files for simulations.
        #
        locs = saH5Py.loadLocalizations("random_map.hdf5")
        locs["z"][:] = 1.0e-3 * settings.z_planes[0]
        saH5Py.saveLocalizations("c1_random_map.hdf5", locs)
        for i in range(1,4):
            locs["x"] += settings.dx
            locs["y"] += settings.dy
            locs["z"][:] = settings.z_planes[i]
            saH5Py.saveLocalizations("c" + str(i+1) + "_random_map.hdf5", locs)

        # Make localization files for simulations.
        #
        locs = saH5Py.loadLocalizations("random_map.hdf5")
        locs["z"][:] = 1.0e-3 * settings.z_planes[0]
        saH5Py.saveLocalizations("c1_random_map.hdf5", locs)
        for i in range(1,4):
            locs["x"] += settings.dx
            locs["y"] += settings.dy
            locs["z"][:] = settings.z_planes[i]
            saH5Py.saveLocalizations("c" + str(i+1) + "_random_map.hdf5", locs)
        
        # Make simulated mapping data.
        # 
        bg_f = lambda s, x, y, h5 : background.UniformBackground(s, x, y, h5, photons = 10)
        cam_f = lambda s, x, y, h5 : camera.SCMOS(s, x, y, h5, "calib.npy")
        pp_f = lambda s, x, y, h5 : photophysics.AlwaysOn(s, x, y, h5, 20000.0)
        psf_f = lambda s, x, y, i3 : psf.GaussianPSF(s, x, y, i3, settings.pixel_size)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)

        for i in range(4):
            sim.simulate("c" + str(i+1) + "_map.dax", "c" + str(i+1) + "_random_map.hdf5", 1)
    
        # Analyze simulated mapping data
        #
        for i in range(4):
            h5_name = "c" + str(i+1) + "_map.hdf5"
            if os.path.exists(h5_name):
                os.remove(h5_name)
            scmos.analyze("c" + str(i+1) + "_map.dax", h5_name, "scmos.xml")
            
        # Measure mapping.
        #
        for i in range(3):
            micrometry.runMicrometry("c1_map.hdf5",
                                     "c" + str(i+2) + "_map.hdf5",
                                     "c1_c" + str(i+2) + "_map.map",
                                     min_size = 5.0,
                                     max_size = 100.0,
                                     max_neighbors = 20,
                                     tolerance = 1.0e-2,
                                     no_plots = True)

        # Merge mapping and save results.
        #
        merged_map = mergeMaps.mergeMaps(["c1_c2_map.map", "c1_c3_map.map", "c1_c4_map.map"])
        
        with open("map.map", 'wb') as fp:
            pickle.dump(merged_map, fp)
        
        # Print mapping.
        #
        if True:
            print("Mapping is:")
            printMapping.printMapping("map.map")
            print("")

        # Check that mapping is close to what we expect (within 5%).
        #
        with open("map.map", 'rb') as fp:
            mappings = pickle.load(fp)

        for i in range(3):
            if not numpy.allclose(mappings["0_" + str(i+1) + "_x"], numpy.array([settings.dx*(i+1), 1.0, 0.0]), rtol = 0.05, atol = 0.05):
                print("X mapping difference for channel", i+1)
            if not numpy.allclose(mappings["0_" + str(i+1) + "_y"], numpy.array([settings.dy*(i+1), 0.0, 1.0]), rtol = 0.05, atol = 0.05):
                print("Y mapping difference for channel", i+1)
    

    ## This part measures / test the PSF measurement.
    ##
    if True:

        # Create drift file, this is used to displace the localizations in the
        # PSF measurement movie.
        #
        dz = numpy.arange(-settings.psf_z_range, settings.psf_z_range + 0.05, 0.01)
        drift_data = numpy.zeros((dz.size, 3))
        drift_data[:,2] = dz
        numpy.savetxt("drift.txt", drift_data)

        # Also create the z-offset file.
        #
        z_offset = numpy.ones((dz.size, 2))
        z_offset[:,1] = dz
        numpy.savetxt("z_offset.txt", z_offset)

        # Create simulated data for PSF measurements.
        #
        bg_f = lambda s, x, y, h5 : background.UniformBackground(s, x, y, h5, photons = 10)
        cam_f = lambda s, x, y, h5 : camera.SCMOS(s, x, y, h5, "calib.npy")
        drift_f = lambda s, x, y, h5 : drift.DriftFromFile(s, x, y, h5, "drift.txt")
        pp_f = lambda s, x, y, h5 : photophysics.AlwaysOn(s, x, y, h5, 20000.0)
        psf_f = lambda s, x, y, h5 : psf.PupilFunction(s, x, y, h5, settings.pixel_size, [])

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                drift_factory = drift_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)

        if True:
            for i in range(4):
                sim.simulate("c" + str(i+1) + "_zcal.dax",
                             "c" + str(i+1) + "_random_map.hdf5",
                             dz.size)
        
        # Get localizations to use for PSF measurement.
        #
        psfLocalizations.psfLocalizations("c1_map_ref.hdf5",
                                          "map.map",
                                          aoi_size = aoi_size)
    
        # Create PSF z stacks.
        #
        for i in range(4):
            psfZStack.psfZStack("c" + str(i+1) + "_zcal.dax",
                                "c1_map_ref_c" + str(i+1) + "_psf.hdf5",
                                "c" + str(i+1) + "_zstack",
                                aoi_size = aoi_size)

        # Measure PSF.
        #
        for i in range(4):
            mpMeasurePSF.measurePSF("c" + str(i+1) + "_zstack.npy",
                                    "z_offset.txt",
                                    "c" + str(i+1) + "_psf_normed.psf",
                                    z_range = settings.psf_z_range,
                                    normalize = True)


    ## This part creates the splines.
    ##
    if True:
        print("Measuring Splines.")
        for i in range(4):
            psfToSpline.psfToSpline("c" + str(i+1) + "_psf_normed.psf",
                                    "c" + str(i+1) + "_psf.spline",
                                    int(settings.psf_size/2))
        
            
    ## This part measures the Cramer-Rao weights.
    ##
    if True:
        print("Calculating weights.")
        planeWeighting.runPlaneWeighting("multicolor.xml",
                                         "weights.npy",
                                         [settings.photons[0][0]],
                                         settings.photons[0][1],
                                         no_plots = True)


if (__name__ == "__main__"):
    configure()
    
