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
import subprocess

import storm_analysis
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.drift as drift
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

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
        params.setAttr("psf0", "filename", "c1_psf_fft.psf")
        params.setAttr("psf1", "filename", "c2_psf_fft.psf")

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
    sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
    subprocess.call(["python", sim_path + "emitters_on_grid.py",
                     "--bin", "grid_list.hdf5",
                     "--nx", str(settings.nx),
                     "--ny", str(settings.ny),
                     "--spacing", "20",
                     "--zrange", str(settings.test_z_range),
                     "--zoffset", str(settings.test_z_offset)])

    # Create randomly located localizations file.
    #
    print("Creating random localization.")
    subprocess.call(["python", sim_path + "emitters_uniform_random.py",
                     "--bin", "random_list.hdf5",
                     "--density", "1.0",
                     "--margin", str(settings.margin),
                     "--sx", str(settings.x_size),
                     "--sy", str(settings.y_size),
                     "--zrange", str(settings.test_z_range)])

    # Create sparser grid for PSF measurement.
    #
    print("Creating data for PSF measurement.")
    sim_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/simulator/"
    subprocess.call(["python", sim_path + "emitters_on_grid.py",
                     "--bin", "psf_list.hdf5",
                     "--nx", "6",
                     "--ny", "3",
                     "--spacing", "40"])

    # Create sCMOS camera calibration files.
    #
    numpy.save("calib.npy", [numpy.zeros((settings.y_size, settings.x_size)) + settings.camera_offset,
                             numpy.ones((settings.y_size, settings.x_size)) * settings.camera_variance,
                             numpy.ones((settings.y_size, settings.x_size)) * settings.camera_gain,
                             1])

    # Create mapping file.
    with open("map.map", 'wb') as fp:
        pickle.dump(settings.mappings, fp)

    if no_splines:
        return

    multiplane_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/multi_plane/"

    # Create pupil functions for 'pupilfn'.
    if (psf_model == "pupilfn"):
        pupilfn_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/pupilfn/"
        print("Creating pupil functions.")
        for i in range(len(settings.z_planes)):
            subprocess.call(["python", pupilfn_path + "make_pupil_fn.py",
                             "--filename", "c" + str(i+1) + "_pupilfn.pfn",
                             "--size", str(settings.psf_size),
                             "--pixel-size", str(settings.pixel_size),
                             "--zmn", str(settings.pupil_fn),
                             "--z-offset", str(-settings.z_planes[i])])

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
        psf_fft_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/psf_fft/"
        spliner_path = os.path.dirname(inspect.getfile(storm_analysis)) + "/spliner/"
        for i in range(len(settings.z_planes)):
            subprocess.call(["python", multiplane_path + "psf_zstack.py",
                             "--movie", "c" + str(i+1) + "_zcal.dax",
                             "--bin", "c" + str(i+1) + "_psf.hdf5",
                             "--zstack", "c" + str(i+1) + "_zstack",
                             "--scmos_cal", "calib.npy",
                             "--aoi_size", str(int(settings.psf_size/2)+1)])

    # Measure PSF and calculate spline for Spliner.
    #
    if (psf_model == "spline"):
    
        # PSFs are independently normalized.
        #
        if settings.independent_heights:
            for i in range(len(settings.z_planes)):
                subprocess.call(["python", multiplane_path + "measure_psf.py",
                                 "--zstack", "c" + str(i+1) + "_zstack.npy",
                                 "--zoffsets", "z_offset.txt",
                                 "--psf_name", "c" + str(i+1) + "_psf_normed.psf",
                                 "--z_range", str(settings.spline_z_range),
                                 "--normalize", "True"])

        # PSFs are normalized to each other.
        #
        else:
            for i in range(len(settings.z_planes)):
                subprocess.call(["python", multiplane_path + "measure_psf.py",
                                 "--zstack", "c" + str(i+1) + "_zstack.npy",
                                 "--zoffsets", "z_offset.txt",
                                 "--psf_name", "c" + str(i+1) + "_psf.psf",
                                 "--z_range", str(settings.spline_z_range)])

            norm_args = ["python", multiplane_path + "normalize_psfs.py",
                         "--psfs", "c1_psf.psf"]
            for i in range(len(settings.z_planes)-1):
                norm_args.append("c" + str(i+2) + "_psf.psf")
            subprocess.call(norm_args)

        # Measure the spline for Spliner.
        #
        print("Measuring Spline.")
        for i in range(len(settings.z_planes)):
            subprocess.call(["python", spliner_path + "psf_to_spline.py",
                             "--psf", "c" + str(i+1) + "_psf_normed.psf",
                             "--spline", "c" + str(i+1) + "_psf.spline",
                             "--spline_size", str(settings.psf_size)])

    # Measure PSF and downsample for PSF FFT.
    #
    elif (psf_model == "psf_fft"):
    
        # PSFs are independently normalized.
        #
        if settings.independent_heights:
            for i in range(len(settings.z_planes)):
                subprocess.call(["python", multiplane_path + "measure_psf.py",
                                 "--zstack", "c" + str(i+1) + "_zstack.npy",
                                 "--zoffsets", "z_offset.txt",
                                 "--psf_name", "c" + str(i+1) + "_psf_normed.psf",
                                 "--z_range", str(settings.psf_z_range),
                                 "--z_step", str(settings.psf_z_step),
                                 "--normalize", "True"])

        # PSFs are normalized to each other.
        #
        else:
            for i in range(len(settings.z_planes)):
                subprocess.call(["python", multiplane_path + "measure_psf.py",
                                 "--zstack", "c" + str(i+1) + "_zstack.npy",
                                 "--zoffsets", "z_offset.txt",
                                 "--psf_name", "c" + str(i+1) + "_psf.psf",
                                 "--z_range", str(settings.psf_z_range),
                                 "--z_step", str(settings.psf_z_step)])

            norm_args = ["python", multiplane_path + "normalize_psfs.py",
                         "--psfs", "c1_psf.psf"]
            for i in range(len(settings.z_planes)-1):
                norm_args.append("c" + str(i+2) + "_psf.psf")
            subprocess.call(norm_args)
        
        # Downsample the PSF to 1x for PSF FFT.
        print("Downsampling PSF.")
        for i in range(len(settings.z_planes)):
            subprocess.call(["python", psf_fft_path + "downsample_psf.py",
                             "--spliner_psf", "c" + str(i+1) + "_psf_normed.psf",
                             "--psf", "c" + str(i+1) + "_psf_fft.psf",
                             "--pixel-size", str(settings.pixel_size)])

    # Calculate Cramer-Rao weighting.
    #
    print("Calculating weights.")
    subprocess.call(["python", multiplane_path + "plane_weighting.py",
                     "--background", str(settings.photons[0][0]),
                     "--photons", str(settings.photons[0][1]),
                     "--output", "weights.npy",
                     "--xml", "multiplane.xml",
                     "--no_plots"])

if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description = 'Multiplane diagnostics configuration.')

    parser.add_argument('--psf-model', dest='psf_model', type=str, required=True,
                        help = "The PSF model, must be one of 'psf_fft', 'pupilfn' or 'spline'")
    parser.add_argument('--no-splines', dest='no_splines', action='store_true', default = False)

    args = parser.parse_args()
    configure(args.psf_model, args.no_splines)
