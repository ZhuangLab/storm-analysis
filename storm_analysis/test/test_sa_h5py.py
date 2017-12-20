#!/usr/bin/env python
"""
Tests of our HDF5 reader, or perhaps more accurately test of our 
understanding of how to use the h5py module.
"""
import numpy
import os
import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5py


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
    with saH5py.SAH5Py(h5_name) as h5:
        h5.addMetadata(metadata)

    # Read metadata.
    with saH5py.SAH5Py(h5_name) as h5:
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
    with saH5py.SAH5Py(h5_name) as h5:
        h5.addLocalizations(peaks, 1)
        h5.addLocalizations(peaks, 1, channel = 1)

    # Read data.
    with saH5py.SAH5Py(h5_name) as h5:

        # Check that frame 0 is empty.
        locs = h5.getLocalizationsInFrame(0)
        assert(not bool(locs))

        # Check frame1.
        locs = h5.getLocalizationsInFrame(1)
        assert(numpy.allclose(peaks["x"], locs["x"]))
        assert(numpy.allclose(peaks["y"], locs["y"]))

        # Check frame1, channel 0.
        locs = h5.getLocalizationsInFrame(1, channel = 1)
        assert(numpy.allclose(peaks["x"], locs["x"]))
        assert(numpy.allclose(peaks["y"], locs["y"]))

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
    with saH5py.SAH5Py(h5_name) as h5:
        h5.addMovieInformation(fr)
        for i in range(fr.getMovieL()):
            h5.addLocalizations(peaks, i)

    # Read data.
    with saH5py.SAH5Py(h5_name) as h5:

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
    with saH5py.SAH5Py(h5_name) as h5:
        h5.addLocalizations(peaks, 1)
        h5.setDriftCorrection(1, dx = 1.0, dy = -1.0)

    # Read data.
    with saH5py.SAH5Py(h5_name) as h5:

        # not corrected.
        locs = h5.getLocalizationsInFrame(1)
        assert(numpy.allclose(peaks["x"], locs["x"]))
        assert(numpy.allclose(peaks["y"], locs["y"]))

        # corrected.
        locs = h5.getLocalizationsInFrame(1, drift_corrected = True)
        assert(numpy.allclose(peaks["x"] + 1.0, locs["x"]))
        assert(numpy.allclose(peaks["y"] - 1.0, locs["y"]))

        
if (__name__ == "__main__"):
    test_sa_h5py_1()
    test_sa_h5py_2()
    test_sa_h5py_3()
    test_sa_h5py_4()


