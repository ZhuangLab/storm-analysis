#!/usr/bin/env python
"""
Make data for testing FRC measurement.

Hazen 01/18
"""
import numpy
import os

import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.diagnostics.frc.settings as settings


def makeData():
    index = 1

    # Create 'tracked' files with different numbers of localizations
    # randomly pulled from the clusters file.
    #
    locs = saH5Py.loadLocalizations("clusters_list.hdf5", fields = ["x", "y"])
    i_arr = numpy.arange(locs["x"].size)

    for reps in settings.n_reps:
    
        wdir = "test_{0:02d}".format(index)
        print(wdir)
        if not os.path.exists(wdir):
            os.makedirs(wdir)
        
        with saH5Py.SAH5Py(os.path.join(wdir, "test.hdf5"), is_existing = False, overwrite = True) as h5:
            h5.setMovieInformation(settings.x_size, settings.y_size, 1, "")
            h5.setPixelSize(settings.pixel_size)
            for i in range(reps):
                track_id = numpy.arange(i * settings.n_tracks_in_group,
                                        (i+1) * settings.n_tracks_in_group,
                                        dtype = numpy.int64)
                numpy.random.shuffle(i_arr)

                lx = locs["x"][i_arr]
                ly = locs["y"][i_arr]
            
                h5.addTracks({"x" : lx[:settings.n_tracks_in_group],
                              "y" : ly[:settings.n_tracks_in_group],
                              "track_id" : track_id})

        index += 1


if (__name__ == "__main__"):
    makeData()
    
            
        
