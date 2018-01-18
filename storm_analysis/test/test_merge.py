#!/usr/bin/env python
"""
Test merging HDF5 files.
"""
import h5py
import numpy
import os

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_utilities.merge_hdf5 as mergeHDF5

    
def test_merge_1():
    """
    Test file merging.
    """
    metadata = "<xml><field1><data1>data</data1></field></xml>"
    ref_tracks = {"x" : numpy.random.randint(0,10,10),
                  "y" : numpy.random.randint(0,10,10)}

    # Create HDF5 files to merge.
    h5_names = []
    for i in range(3):
        h5_name = storm_analysis.getPathOutputTest("test_merge_f" + str(i) + ".hdf5")
        h5_names.append(h5_name)

        with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
            h5.setMovieInformation(20,20,1,"")
            h5.addMetadata(metadata)
            h5.addTracks(ref_tracks)
            
    # Merge.
    merge_name = storm_analysis.getPathOutputTest("test_merge.hdf5")
    storm_analysis.removeFile(merge_name)
    mergeHDF5.mergeHDF5(h5_names, merge_name)

    # Check merge.
    with saH5Py.SAH5Py(merge_name) as h5:
        assert(metadata == h5.getMetadata())
        for tracks in h5.tracksIterator():
            assert(numpy.allclose(ref_tracks["x"], tracks["x"]))
            

def test_merge_2():
    """
    Test file merging, skipping files with no tracks.
    """
    metadata = "<xml><field1><data1>data</data1></field></xml>"
    ref_tracks = {"x" : numpy.random.randint(0,10,10),
                  "y" : numpy.random.randint(0,10,10)}

    # Create HDF5 files to merge.
    h5_names = []
    for i in range(3):
        h5_name = storm_analysis.getPathOutputTest("test_merge_f" + str(i) + ".hdf5")
        h5_names.append(h5_name)

        with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
            h5.setMovieInformation(20,20,1,"")
            h5.addMetadata(metadata)

            if(i != 1):
                h5.addTracks(ref_tracks)

    # Merge.
    merge_name = storm_analysis.getPathOutputTest("test_merge.hdf5")
    storm_analysis.removeFile(merge_name)
    mergeHDF5.mergeHDF5(h5_names, merge_name)

    # Check merge.
    with saH5Py.SAH5Py(merge_name) as h5:
        assert(metadata == h5.getMetadata())
        for tracks in h5.tracksIterator():
            assert(numpy.allclose(ref_tracks["x"], tracks["x"]))            


if (__name__ == "__main__"):
    test_merge_1()
    test_merge_2()
