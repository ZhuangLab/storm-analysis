#!/usr/bin/env python
"""
Make data for testing SLURM scripts.

Hazen 09/18
"""
import glob
import numpy
import os
import pickle

import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

import storm_analysis.slurm.split_analysis_xml as splitAnalysisXML

import storm_analysis.diagnostics.slurm.settings as settings


def makeData():

    # Create .bin files for each plane.
    h5_locs = saH5Py.loadLocalizations("grid_list.hdf5")
        
    # Load channel to channel mapping file.
    with open("map.map", 'rb') as fp:
        mappings = pickle.load(fp)

    # Add z offset to reference localizations.
    x = h5_locs["x"].copy()
    y = h5_locs["y"].copy()
    z = h5_locs["z"].copy() + settings.z_planes[0]
    h5_temp = {"x" : x,
               "y" : y,
               "z" : z}
    saH5Py.saveLocalizations("sim_input_c1.hdf5", h5_temp)


    # Create a movie for first plane.
    [bg, photons] = settings.photons

    # Adjust photons by the number of planes.
    photons = photons/float(len(settings.z_planes))

    bg_f = lambda s, x, y, i3 : background.UniformBackground(s, x, y, i3, photons = bg)
    cam_f = lambda s, x, y, i3 : camera.SCMOS(s, x, y, i3, "calib.npy")
    pp_f = lambda s, x, y, i3 : photophysics.SimpleSTORM(s, x, y, i3,
                                                         photons = photons,
                                                         on_time = settings.on_time,
                                                         off_time = settings.off_time)
    psf_f = lambda s, x, y, i3 : psf.PupilFunction(s, x, y, i3, settings.pixel_size, settings.pupil_fn)

    sim = simulate.Simulate(background_factory = bg_f,
                            camera_factory = cam_f,
                            photophysics_factory = pp_f,
                            psf_factory = psf_f,
                            x_size = settings.x_size,
                            y_size = settings.y_size)

    sim.simulate(os.path.join(settings.wdir, "test_c1.dax"),
                 "sim_input_c1.hdf5",
                 settings.n_frames)

    # Create other movies.
    for i in range(1, len(settings.z_planes)):
        cx = mappings["0_" + str(i) + "_x"]
        cy = mappings["0_" + str(i) + "_y"]
        z_offset = settings.z_planes[i] - settings.z_planes[0]
            
        pp_f = lambda s, x, y, i3 : photophysics.Duplicate(s, x, y, i3,
                                                           h5_name = os.path.join(settings.wdir, "test_c1_ref.hdf5"),
                                                           cx = cx, cy = cy, z_offset = z_offset)

        sim = simulate.Simulate(background_factory = bg_f,
                                camera_factory = cam_f,
                                photophysics_factory = pp_f,
                                psf_factory = psf_f,
                                x_size = settings.x_size,
                                y_size = settings.y_size)

        sim.simulate(os.path.join(settings.wdir, "test_c" + str(i+1) + ".dax"),
                     "sim_input_c1.hdf5", # This is not actually used.
                     settings.n_frames)            

    # Remove any old XML files.
    for elt in glob.glob(os.path.join(settings.wdir, "job*.xml")):
        os.remove(elt)
            
    # Make analysis XML files.
    splitAnalysisXML.splitAnalysisXML(settings.wdir, "multiplane.xml", 0, settings.n_frames, settings.divisions)
        

if (__name__ == "__main__"):
    makeData()
    
