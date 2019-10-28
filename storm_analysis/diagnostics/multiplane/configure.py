#!/usr/bin/env python
"""
Configure folder for Multiplane testing.

Hazen 10/17
"""
import argparse
import inspect
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

import storm_analysis.multi_plane.measure_psf as mpMeasurePSF
import storm_analysis.multi_plane.normalize_psfs as normalizePSFs
import storm_analysis.multi_plane.plane_weighting as planeWeighting
import storm_analysis.multi_plane.psf_zstack as psfZStack

import storm_analysis.pupilfn.make_pupil_fn as makePupilFn

import storm_analysis.spliner.psf_to_spline as psfToSpline

import storm_analysis.diagnostics.multiplane.settings as settings


def testingParameters(psf_model):
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

    if (psf_model == "psf_fft"):
        params.setAttr("psf0", "filename", "c1_psf_normed.psf")
        params.setAttr("psf1", "filename", "c2_psf_normed.psf")

    elif (psf_model == "pupilfn"):
        params.setAttr("pupilfn0", "filename", "c1_pupilfn.pfn")
        params.setAttr("pupilfn1", "filename", "c2_pupilfn.pfn")
        
    elif (psf_model == "spline"):
        params.setAttr("spline0", "filename", "c1_psf.spline")
        params.setAttr("spline1", "filename", "c2_psf.spline")
        
    else:
        raise Exception("Unknown PSF model " + args.psf_model)

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.0")
    
    if (psf_model == "psf_fft"):
        params.setAttr("max_z", "float", str(settings.psf_z_range + 0.001))
        params.setAttr("min_z", "float", str(-(settings.psf_z_range - 0.001)))

    elif (psf_model == "pupilfn"):
        params.setAttr("max_z", "float", str(settings.pupilfn_z_range))
        params.setAttr("min_z", "float", str(-settings.pupilfn_z_range))
        
    elif (psf_model == "spline"):
        params.setAttr("max_z", "float", str(settings.spline_z_range + 0.001))
        params.setAttr("min_z", "float", str(-(settings.spline_z_range - 0.001)))

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    # 'peak_locations' testing.
    if hasattr(settings, "peak_locations") and (settings.peak_locations is not None):
        params.setAttr("peak_locations", "filename", settings.peak_locations)    

    return params

def configure(psf_model, no_splines):
    # Create parameters file for analysis.
    #
    print("Creating XML file.")
    params = testingParameters(psf_model)
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

    # Create randomly located localizations file.
    #
    print("Creating random localization.")
    emittersUniformRandom.emittersUniformRandom("random_list.hdf5",
                                                1.0,
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

    # Create sCMOS camera calibration files.
    #
    numpy.save("calib.npy", [numpy.zeros((settings.y_size, settings.x_size)) + settings.camera_offset,
                             numpy.ones((settings.y_size, settings.x_size)) * settings.camera_variance,
                             numpy.ones((settings.y_size, settings.x_size)) * settings.camera_gain,
                             numpy.ones((settings.y_size, settings.x_size)),
                             2])

    # Create mapping file.
    with open("map.map", 'wb') as fp:
        pickle.dump(settings.mappings, fp)

    if no_splines:
        return

    # Create pupil functions for 'pupilfn'.
    if (psf_model == "pupilfn"):
        print("Creating pupil functions.")
        for i in range(len(settings.z_planes)):
            makePupilFn.makePupilFunction("c" + str(i+1) + "_pupilfn.pfn",
                                          settings.psf_size,
                                          settings.pixel_size * 1.0e-3,
                                          settings.pupil_fn,
                                          z_offset = -settings.z_planes[i])


    # Both 'spline' and 'psf_fft' need measured PSFs.
    else:
    
        # Create localization files for PSF measurement.
        #
        locs = saH5Py.loadLocalizations("psf_list.hdf5")

        for i, z_offset in enumerate(settings.z_planes):
            cx = settings.mappings["0_" + str(i) + "_x"]
            cy = settings.mappings["0_" + str(i) + "_y"]
            locs_temp = {"x" : locs["x"].copy(),
                         "y" : locs["y"].copy(),
                         "z" : locs["z"].copy()}
            xi = locs_temp["x"]
            yi = locs_temp["y"]
            xf = cx[0] + cx[1] * xi + cx[2] * yi
            yf = cy[0] + cy[1] * xi + cy[2] * yi
            locs_temp["x"] = xf
            locs_temp["y"] = yf
            locs_temp["z"][:] = z_offset

            saH5Py.saveLocalizations("c" + str(i+1) + "_psf.hdf5", locs_temp)

        # Create drift file, this is used to displace the localizations in the
        # PSF measurement movie.
        #
        dz = numpy.arange(-settings.spline_z_range, settings.spline_z_range + 0.001, 0.01)
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
        psf_f = lambda s, x, y, h5 : psf.PupilFunction(s, x, y, h5, settings.pixel_size, settings.pupil_fn)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                drift_factory = drift_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)

        for i in range(len(settings.z_planes)):
            sim.simulate("c" + str(i+1) + "_zcal.dax",
                         "c" + str(i+1) + "_psf.hdf5",
                         dz.size)
        
        # Measure the PSF.
        #
        print("Measuring PSFs.")
        for i in range(len(settings.z_planes)):
            psfZStack.psfZStack("c" + str(i+1) + "_zcal.dax",
                                "c" + str(i+1) + "_psf.hdf5",
                                "c" + str(i+1) + "_zstack",
                                aoi_size = int(settings.psf_size/2 + 1))
            

    # Measure PSF and calculate spline for Spliner.
    #
    if (psf_model == "spline"):
    
        # PSFs are independently normalized.
        #
        if settings.independent_heights:
            for i in range(len(settings.z_planes)):
                mpMeasurePSF.measurePSF("c" + str(i+1) + "_zstack.npy",
                                        "z_offset.txt",
                                        "c" + str(i+1) + "_psf_normed.psf",
                                        z_range = settings.spline_z_range,
                                        normalize = True)

        # PSFs are normalized to each other.
        #
        else:
            for i in range(len(settings.z_planes)):
                mpMeasurePSF.measurePSF("c" + str(i+1) + "_zstack.npy",
                                        "z_offset.txt",
                                        "c" + str(i+1) + "_psf.psf",
                                        z_range = settings.spline_z_range)


            norm_args = ["c1_psf.psf"]
            for i in range(len(settings.z_planes)-1):
                norm_args.append("c" + str(i+2) + "_psf.psf")
            normalizePSFs.normalizePSFs(norm_args)

        # Measure the spline for Spliner.
        #
        print("Measuring Spline.")
        for i in range(len(settings.z_planes)):
            psfToSpline.psfToSpline("c" + str(i+1) + "_psf_normed.psf",
                                    "c" + str(i+1) + "_psf.spline",
                                    int(settings.psf_size/2))


    # Measure PSF and downsample for PSF FFT.
    #
    elif (psf_model == "psf_fft"):
    
        # PSFs are independently normalized.
        #
        if settings.independent_heights:
            for i in range(len(settings.z_planes)):
                mpMeasurePSF.measurePSF("c" + str(i+1) + "_zstack.npy",
                                        "z_offset.txt",
                                        "c" + str(i+1) + "_psf_normed.psf",
                                        z_range = settings.spline_z_range,
                                        normalize = True)

        # PSFs are normalized to each other.
        #
        else:
            for i in range(len(settings.z_planes)):
                mpMeasurePSF.measurePSF("c" + str(i+1) + "_zstack.npy",
                                        "z_offset.txt",
                                        "c" + str(i+1) + "_psf.psf",
                                        z_range = settings.spline_z_range)


            norm_args = ["c1_psf.psf"]
            for i in range(len(settings.z_planes)-1):
                norm_args.append("c" + str(i+2) + "_psf.psf")
            normalizePSFs.normalizePSFs(norm_args)


    # Calculate Cramer-Rao weighting.
    #
    print("Calculating weights.")
    planeWeighting.runPlaneWeighting("multiplane.xml",
                                     "weights.npy",
                                     [settings.photons[0][0]],
                                     settings.photons[0][1],
                                     no_plots = True)


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description = 'Multiplane diagnostics configuration.')

    parser.add_argument('--psf-model', dest='psf_model', type=str, required=True,
                        help = "The PSF model, must be one of 'psf_fft', 'pupilfn' or 'spline'")
    parser.add_argument('--no-splines', dest='no_splines', action='store_true', default = False)

    args = parser.parse_args()
    configure(args.psf_model, args.no_splines)
