#!/usr/bin/env python
"""
Test merging HDF5 files.
"""
import h5py
import numpy
import os

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_utilities.align_and_merge as alignAndMerge
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
            h5.setPixelSize(100.0)
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
            h5.addMetadata(metadata)
            h5.setMovieInformation(20,20,1,"")
            h5.setPixelSize(100.0)
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


def test_align_merge_1():
    """
    Test aligning and merging two HDF5 files.
    """
    n_locs = 500
    tracks = {"x" : numpy.random.normal(loc = 10.0, scale = 0.2, size = n_locs),
              "y" : numpy.random.normal(loc = 10.0, scale = 0.2, size = n_locs),
              "z" : numpy.random.normal(scale = 0.05, size = n_locs)}

    h5_in1 = storm_analysis.getPathOutputTest("test_align_merge_1.hdf5")
    h5_in2 = storm_analysis.getPathOutputTest("test_align_merge_2.hdf5")
    h5_alm = storm_analysis.getPathOutputTest("test_align_merge_3.hdf5")

    # Create input files.
    t_dx = 1.0
    t_dz = 0.3
    with saH5Py.SAH5Py(h5_in1, is_existing = False, overwrite = True) as h5:
        h5.addMetadata("<xml><field1><data1>1</data1></field></xml>")
        h5.setMovieInformation(20, 20, 2, "")
        h5.setPixelSize(100.0)
        h5.addTracks(tracks)

    with saH5Py.SAH5Py(h5_in2, is_existing = False, overwrite = True) as h5:
        h5.addMetadata("<xml><field1><data1>2</data1></field></xml>")
        h5.setMovieInformation(20, 20, 2, "")
        h5.setPixelSize(100.0)

        tracks["x"] += t_dx
        tracks["z"] += t_dz
        h5.addTracks(tracks)

    # Align and merge.
    storm_analysis.removeFile(h5_alm)
    [dx, dy, dz] = alignAndMerge.alignAndMerge(h5_in1, h5_in2, h5_alm)

    # Check that we got the right offsets.
    assert(numpy.allclose(numpy.array([dx, dy, dz]),
                          numpy.array([-t_dx, 0.0, -t_dz]),
                          atol = 0.001,
                          rtol = 0.1))

    # Check that the output file is correctly aligned.
    with saH5Py.SAH5Py(h5_alm) as h5:
        tracks = h5.getTracks(fields = ["x", "y", "z"])
        assert(numpy.allclose(numpy.array([numpy.std(tracks["x"]),
                                           numpy.std(tracks["y"]),
                                           numpy.std(tracks["z"])]),
                              numpy.array([0.2, 0.2, 0.05]),
                              atol = 0.001,
                              rtol = 0.1))

    
def test_align_merge_2():
    """
    Test aligning and merging two HDF5 files with offset.
    """
    n_locs = 500
    tracks = {"x" : numpy.random.normal(loc = 10.0, scale = 0.2, size = n_locs),
              "y" : numpy.random.normal(loc = 10.0, scale = 0.2, size = n_locs),
              "z" : numpy.random.normal(scale = 0.05, size = n_locs)}

    h5_in1 = storm_analysis.getPathOutputTest("test_align_merge_1.hdf5")
    h5_in2 = storm_analysis.getPathOutputTest("test_align_merge_2.hdf5")
    h5_alm = storm_analysis.getPathOutputTest("test_align_merge_3.hdf5")

    # Create input files.
    t_dx = 2.0
    t_dz = 0.3
    with saH5Py.SAH5Py(h5_in1, is_existing = False, overwrite = True) as h5:
        h5.addMetadata("<xml><field1><data1>1</data1></field></xml>")
        h5.setMovieInformation(20, 20, 2, "")
        h5.setPixelSize(100.0)
        h5.addTracks(tracks)

    with saH5Py.SAH5Py(h5_in2, is_existing = False, overwrite = True) as h5:
        h5.addMetadata("<xml><field1><data1>2</data1></field></xml>")
        h5.setMovieInformation(20, 20, 2, "")
        h5.setPixelSize(100.0)

        tracks["x"] += t_dx
        tracks["z"] += t_dz
        h5.addTracks(tracks)

    # Align and merge with offset.
    storm_analysis.removeFile(h5_alm)
    [dx, dy, dz] = alignAndMerge.alignAndMerge(h5_in1, h5_in2, h5_alm, dx = -t_dx)

    # Check that we got the right offsets.
    assert(numpy.allclose(numpy.array([dx, dy, dz]),
                          numpy.array([-t_dx, 0.0, -t_dz]),
                          atol = 0.001,
                          rtol = 0.1))


    # Check that the output file is correctly aligned.
    with saH5Py.SAH5Py(h5_alm) as h5:
        tracks = h5.getTracks(fields = ["x", "y", "z"])
        assert(numpy.allclose(numpy.array([numpy.std(tracks["x"]),
                                           numpy.std(tracks["y"]),
                                           numpy.std(tracks["z"])]),
                              numpy.array([0.2, 0.2, 0.05]),
                              atol = 0.001,
                              rtol = 0.1))
        

if (__name__ == "__main__"):
    test_merge_1()
    test_merge_2()
    test_align_merge_1()
    test_align_merge_2()
