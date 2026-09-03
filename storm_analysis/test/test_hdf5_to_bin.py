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

    
def test_hdf5_to_bin_track_normalization():
    """
    Tracks store most of their fields as sums over the localizations in the
    track, so exporting them raw made a molecule that stayed on for ten
    frames ten times as wide.

    Width, aspect ratio and background must not depend on the track length.
    'a' and 'h' must, they are the total photons and the total height over
    the track, which is how the Insight3 era averager reported them.
    """
    [sx, sy] = [1.5, 2.0]
    [height, background, photons] = [1000.0, 20.0, 5000.0]
    pixel_size = 100.0

    lengths = numpy.array([1, 3, 10], dtype = numpy.int32)
    n = lengths.size

    tracks = {"x" : numpy.arange(n, dtype = numpy.float64) + 10.0,
              "y" : numpy.arange(n, dtype = numpy.float64) + 10.0,
              "z" : numpy.zeros(n),
              "category" : numpy.zeros(n, dtype = numpy.int32),
              "frame_number" : numpy.ones(n, dtype = numpy.int32),
              "track_id" : numpy.arange(n, dtype = numpy.int64),
              "track_length" : lengths,
              "xsigma" : lengths*sx,
              "ysigma" : lengths*sy,
              "height" : lengths*height,
              "background" : lengths*background,
              "sum" : lengths*photons}

    h5_name = storm_analysis.getPathOutputTest("test_h5_to_bin_norm.hdf5")
    storm_analysis.removeFile(h5_name)

    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addMetadata("<settings/>")
        h5.setMovieInformation(256, 256, 10, "XYZZY")
        h5.setPixelSize(pixel_size)
        h5.addLocalizations({"x" : numpy.zeros(1), "y" : numpy.zeros(1)}, 0)
        h5.addTracks(tracks)

    i3_name = storm_analysis.getPathOutputTest("test_h5_to_bin_norm.bin")
    storm_analysis.removeFile(i3_name)
    hdf5ToBin.hdf5ToBin(h5_name, i3_name)

    i3_data = readinsight3.loadI3File(i3_name, verbose = False)

    # Insight3 widths are 2 * sigma, in nanometers.
    [wx, wy] = [2.0*sx*pixel_size, 2.0*sy*pixel_size]

    assert(numpy.allclose(i3_data['w'], numpy.sqrt(wx*wy)*numpy.ones(n)))
    assert(numpy.allclose(i3_data['ax'], (wy/wx)*numpy.ones(n)))
    assert(numpy.allclose(i3_data['bg'], background*numpy.ones(n)))

    # Totals, so these two do scale with the track length.
    assert(numpy.allclose(i3_data['a'], lengths*photons))
    assert(numpy.allclose(i3_data['h'], lengths*height))


if (__name__ == "__main__"):
    test_hdf5_to_bin_1()
    test_hdf5_to_bin_2()
    test_hdf5_to_bin_track_normalization()

