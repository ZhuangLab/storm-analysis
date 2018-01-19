#!/usr/bin/env python
"""
Tests of the HDF5 to bin converter, which is probably going to be
important, at least for a while.
"""
import numpy

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_utilities.hdf5_to_bin as hdf5ToBin


def test_hdf5_to_bin_1():
    """
    Test localizations conversion.
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    h5_name = storm_analysis.getPathOutputTest("test_sa_hdf5.hdf5")
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addMetadata("<settings/>")
        h5.setMovieInformation(256, 256, 10, "XYZZY")
        h5.setPixelSize(100.0)
        h5.addLocalizations(peaks, 1)

    # Convert.
    i3_name = storm_analysis.getPathOutputTest("test_mlist.bin")
    storm_analysis.removeFile(i3_name)
    hdf5ToBin.hdf5ToBin(h5_name, i3_name)

    # Load Insight3 file and check values.
    i3_data = readinsight3.loadI3File(i3_name, verbose = False)

    assert(numpy.allclose(peaks["x"], i3_data['x'] - 1.0))
    assert(numpy.allclose(peaks["y"], i3_data['y'] - 1.0))
    assert(numpy.allclose(i3_data['fr'], 2*numpy.ones(10)))


def test_hdf5_to_bin_2():
    """
    Test tracks conversion.
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    h5_name = storm_analysis.getPathOutputTest("test_sa_hdf5.hdf5")
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addMetadata("<settings/>")
        h5.setMovieInformation(256, 256, 10, "XYZZY")
        h5.setPixelSize(100.0)
        h5.addTracks(peaks)

    # Convert.
    i3_name = storm_analysis.getPathOutputTest("test_mlist.bin")
    storm_analysis.removeFile(i3_name)
    hdf5ToBin.hdf5ToBin(h5_name, i3_name)

    # Load Insight3 file and check values.
    i3_data = readinsight3.loadI3File(i3_name, verbose = False)

    assert(numpy.allclose(peaks["x"], i3_data['x'] - 1.0))
    assert(numpy.allclose(peaks["y"], i3_data['y'] - 1.0))    
    assert(numpy.allclose(i3_data['fr'], numpy.ones(10)))

    
if (__name__ == "__main__"):
    test_hdf5_to_bin_1()
    test_hdf5_to_bin_2()

