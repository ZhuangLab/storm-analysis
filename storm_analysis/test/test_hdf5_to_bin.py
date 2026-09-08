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


def test_hdf5_to_bin_track_expansion():
    """
    Tracks carry their start frame and length through to the Insight3 file,
    and bin_to_lmchallenge_format expands them back into one row per frame.

    hdf5ToBin() wrote every track into frame 1 and left 'tl' at its default
    of 1, so the expansion was silently a no-op and a movie's worth of
    tracks came out as one row each, all in frame 1.
    """
    import storm_analysis.sa_utilities.bin_to_lmchallenge_format as binToLMC

    photons = 6000.0
    pixel_size = 100.0

    lengths = numpy.array([1, 3, 10], dtype = numpy.int32)
    starts = numpy.array([0, 4, 20], dtype = numpy.int32)
    n = lengths.size

    tracks = {"x" : numpy.arange(n, dtype = numpy.float64) + 10.0,
              "y" : numpy.arange(n, dtype = numpy.float64) + 20.0,
              "z" : numpy.zeros(n),
              "category" : numpy.zeros(n, dtype = numpy.int32),
              "frame_number" : starts,
              "track_id" : numpy.arange(n, dtype = numpy.int64),
              "track_length" : lengths,
              "xsigma" : lengths*1.5,
              "ysigma" : lengths*1.5,
              "height" : lengths*1000.0,
              "background" : lengths*20.0,
              "sum" : lengths*photons}

    h5_name = storm_analysis.getPathOutputTest("test_h5_to_bin_expand.hdf5")
    storm_analysis.removeFile(h5_name)

    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addMetadata("<settings/>")
        h5.setMovieInformation(256, 256, 40, "XYZZY")
        h5.setPixelSize(pixel_size)
        h5.addLocalizations({"x" : numpy.zeros(1), "y" : numpy.zeros(1)}, 0)
        h5.addTracks(tracks)

    i3_name = storm_analysis.getPathOutputTest("test_h5_to_bin_expand.bin")
    storm_analysis.removeFile(i3_name)
    hdf5ToBin.hdf5ToBin(h5_name, i3_name)

    i3_data = readinsight3.loadI3File(i3_name, verbose = False)

    # The track length survives, and the frame is where the track started.
    # Insight3 frame numbers start at 1.
    assert(numpy.allclose(i3_data['tl'], lengths))
    assert(numpy.allclose(i3_data['fr'], starts + 1))

    # 'a' is still the whole track's photons.
    assert(numpy.allclose(i3_data['a'], lengths*photons))

    ## Now expand it.
    txt_name = storm_analysis.getPathOutputTest("test_h5_to_bin_expand.txt")
    storm_analysis.removeFile(txt_name)

    n_rows = binToLMC.binToLMChallenge(i3_name, txt_name, pixel_size, verbose = False)

    assert(n_rows == numpy.sum(lengths))

    with open(txt_name) as fp:
        rows = [line.split(",") for line in fp.readlines()[1:]]

    assert(len(rows) == numpy.sum(lengths))

    # Each track occupies consecutive frames from where it started, and each
    # row carries its own share of the photons rather than the whole track's.
    at = 0
    for i in range(n):
        for j in range(lengths[i]):
            frame = int(rows[at][1])
            intensity = float(rows[at][5])
            assert(frame == starts[i] + 1 + j)
            assert(abs(intensity - photons) < 1.0e-3)
            at += 1


if (__name__ == "__main__"):
    test_hdf5_to_bin_1()
    test_hdf5_to_bin_2()
    test_hdf5_to_bin_track_normalization()
    test_hdf5_to_bin_track_expansion()

