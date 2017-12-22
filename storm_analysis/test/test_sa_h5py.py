#!/usr/bin/env python
"""
Tests of our HDF5 reader, or perhaps more accurately test of our 
understanding of how to use the h5py module.
"""
import h5py
import numpy
import os

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_library.writeinsight3 as i3w


class FakeReader(object):
    """
    Fake movie reader object.
    """
    def __init__(self, **kwds):
        super(FakeReader, self).__init__(**kwds)

    def hashID(self):
        return "xyzzy"

    def getMovieL(self):
        return 10

    def getMovieX(self):
        return 128

    def getMovieY(self):
        return 128

    
def test_sa_h5py_1():
    """
    Test metadata round trip.
    """
    metadata = "<xml><field1><data1>data</data1></field></xml>"

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write metadata.
    with saH5Py.SAH5Py(h5_name) as h5:
        h5.addMetadata(metadata)

    # Read metadata.
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(metadata == h5.getMetadata())


def test_sa_h5py_2():
    """
    Test data round trip.
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name) as h5:
        h5.addLocalizations(peaks, 1)
        h5.addLocalizations(peaks, 1, channel = 1)

    # Read data.
    with saH5Py.SAH5Py(h5_name) as h5:

        # Check that frame 0 is empty.
        locs = h5.getLocalizationsInFrame(0)
        assert(not bool(locs))

        # Check frame1.
        locs = h5.getLocalizationsInFrame(1)
        assert(numpy.allclose(peaks["x"], locs["x"]))
        assert(numpy.allclose(peaks["y"], locs["y"]))
        assert(numpy.allclose(peaks["x"], locs["c1_x"]))
        assert(numpy.allclose(peaks["y"], locs["c1_y"]))

        # Check getting a specific field.
        locs = h5.getLocalizationsInFrame(1, fields = ["x"])
        assert("x" in locs)
        assert(not "y" in locs)


def test_sa_h5py_3():
    """
    Test getting data from multiple frames.
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    fr = FakeReader()
    with saH5Py.SAH5Py(h5_name) as h5:
        h5.addMovieInformation(fr)
        for i in range(fr.getMovieL()):
            h5.addLocalizations(peaks, i)

    # Read data.
    with saH5Py.SAH5Py(h5_name) as h5:

        # Check localizations in first 5 frames.
        locs = h5.getLocalizationsInFrameRange(0,5)
        assert(locs["x"].size == 50)

        # Get all the localizations.
        locs = h5.getLocalizations()
        assert(locs["x"].size == (10.0 * fr.getMovieL()))
        

def test_sa_h5py_4():
    """
    Test handling of drift correction.
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name) as h5:
        h5.addLocalizations(peaks, 1)
        h5.setDriftCorrection(1, dx = 1.0, dy = -1.0)

    # Read data.
    with saH5Py.SAH5Py(h5_name) as h5:

        # not corrected.
        locs = h5.getLocalizationsInFrame(1)
        assert(numpy.allclose(peaks["x"], locs["x"]))
        assert(numpy.allclose(peaks["y"], locs["y"]))

        # corrected.
        locs = h5.getLocalizationsInFrame(1, drift_corrected = True)
        assert(numpy.allclose(peaks["x"] + 1.0, locs["x"]))
        assert(numpy.allclose(peaks["y"] - 1.0, locs["y"]))


def test_sa_h5py_5():
    """
    Test querying if the HDF5 file is a storm-analysis file.
    """
    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Open empty file.
    with saH5Py.SAH5Py(h5_name) as h5:
        pass
    assert(saH5Py.isSAHDF5(h5_name))

    # Create generic HDF5 file.
    f = h5py.File(h5_name, "w")
    f.close()
    assert(not saH5Py.isSAHDF5(h5_name))

    # Create Insight3 file.
    with i3w.I3Writer(h5_name) as i3:
        pass
    assert not(saH5Py.isSAHDF5(h5_name))


def test_sa_h5py_6():
    """
    Test adding tracks.
    """
    tracks = {"x" : numpy.zeros(10),
              "y" : numpy.ones(10)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write tracks.
    with saH5Py.SAH5Py(h5_name) as h5:
        h5.addTracks(tracks)

    # Read tracks.
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.getNTracks() == tracks["x"].size)

    # Write tracks again, this should overwrite above.
    with saH5Py.SAH5Py(h5_name) as h5:
        h5.addTracks(tracks)
        h5.addTracks(tracks)

    # Read tracks.
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.getNTracks() == 2*tracks["x"].size)


def test_sa_h5py_7():
    """
    Test tracks iterator.
    """
    tracks = {"x" : numpy.zeros(10),
              "y" : numpy.ones(10)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # No tracks.
    with saH5Py.SAH5Py(h5_name) as h5:
        pass

    with saH5Py.SAH5Py(h5_name) as h5:
        for t in h5.tracksIterator():
            assert(False) # We should not get here.

    # Tracks.
    with saH5Py.SAH5Py(h5_name) as h5:
        h5.addTracks(tracks)

    with saH5Py.SAH5Py(h5_name) as h5:
        for t in h5.tracksIterator():
            assert(numpy.allclose(t["x"], tracks["x"]))

        # Only get one field.
        for t in h5.tracksIterator(["x"]):
            assert(not "y" in t)



def test_sa_h5py_8():
    """
    Test localizations iterator.
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    fr = FakeReader()
    with saH5Py.SAH5Py(h5_name) as h5:
        h5.addMovieInformation(fr)
        for i in range(fr.getMovieL()):
            h5.addLocalizations(peaks, 2*i)

    # Read data.
    with saH5Py.SAH5Py(h5_name) as h5:
        for locs in h5.localizationsIterator():
            assert(locs["x"].size == 10)


if (__name__ == "__main__"):
    test_sa_h5py_1()
    test_sa_h5py_2()
    test_sa_h5py_3()
    test_sa_h5py_4()
    test_sa_h5py_5()
    test_sa_h5py_6()
    test_sa_h5py_7()
    test_sa_h5py_8()
