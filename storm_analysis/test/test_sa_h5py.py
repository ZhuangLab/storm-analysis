#!/usr/bin/env python
"""
Tests of our HDF5 reader, or perhaps more accurately test of our 
understanding of how to use the h5py module.
"""
import h5py
import numpy

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
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
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
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.setMovieInformation(1,1,2,"")
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
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addMovieInformation(fr)
        for i in range(fr.getMovieL()):
            h5.addLocalizations(peaks, i)

    # Read data.
    with saH5Py.SAH5Reader(h5_name) as h5:

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
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
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
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
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
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
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
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        pass

    with saH5Py.SAH5Py(h5_name) as h5:
        for t in h5.tracksIterator():
            assert(False) # We should not get here.

    # Tracks.
    storm_analysis.removeFile(h5_name)
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
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
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addMovieInformation(fr)
        for i in range(fr.getMovieL()):
            h5.addLocalizations(peaks, 2*i)

    # Read data.
    with saH5Py.SAH5Py(h5_name) as h5:
        for fnum, locs in h5.localizationsIterator():
            assert(locs["x"].size == 10)


def test_sa_h5py_9():
    """
    Test setting the track id field.
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Add localizations and track id.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addLocalizations(peaks, 1)
        h5.addTrackID(numpy.ones(10), 1)

    # Check track id.
    with saH5Py.SAH5Py(h5_name) as h5:
        locs = h5.getLocalizationsInFrame(1)
        assert(numpy.allclose(locs["track_id"], numpy.ones(10)))
        
    # Change track id.
    with saH5Py.SAH5Py(h5_name) as h5:
        h5.addTrackID(numpy.zeros(10), 1)

    # Check track id.
    with saH5Py.SAH5Py(h5_name) as h5:
        locs = h5.getLocalizationsInFrame(1)
        assert(numpy.allclose(locs["track_id"], numpy.zeros(10)))


def test_sa_h5py_10():
    """
    Test 'is_existing' and 'overwrite' parameters.
    """
    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Test failure on trying to open a file that does not exist.
    try:
        with saH5Py.SAH5Py(h5_name) as h5:
            pass
    except saH5Py.SAH5PyException:
        pass
    else:
        assert(False)

    # Test failure on trying to overwrite a file that does exist.

    # Create the file.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        pass

    # Test that we cannot overwrite it.
    try:
        with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
            pass
    except saH5Py.SAH5PyException:
        pass
    else:
        assert(False)

    # Test that we can overwrite it.
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        pass


def test_sa_h5py_11():
    """
    Test hasLocalizationField() and hasTracksField()
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.setMovieInformation(256, 256, 10, "XYZZY")
        h5.addLocalizations(peaks, 1)
        h5.addTracks(peaks)

    # Check.
    with saH5Py.SAH5Py(h5_name) as h5:

        assert(h5.hasLocalizationsField("x"))
        assert(not h5.hasLocalizationsField("x1"))

        assert(h5.hasTracksField("x"))
        assert(not h5.hasTracksField("x1"))


def test_sa_h5py_12():
    """
    Test handling of multiple channels.
    """
    peaks = {"x" : numpy.zeros(3),
             "y" : numpy.ones(3)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.setMovieInformation(1,1,2,"")
        h5.addLocalizations(peaks, 1)

        peaks["x"] += 1
        peaks["y"] += 1
        h5.addLocalizations(peaks, 1, channel = 1)

        peaks["x"] += 1
        peaks["y"] += 1        
        h5.addLocalizations(peaks, 1, channel = 2)

    # Read data.
    with saH5Py.SAH5Py(h5_name) as h5:

        # Check getting number of channels.
        assert(h5.getNChannels() == 3)

        for [fnum, locs] in h5.localizationsIterator():
            for i, elt in enumerate(h5.splitByChannel(locs)):
                assert(numpy.allclose(elt["x"], i * numpy.ones(3)))
                assert(numpy.allclose(elt["y"], i * numpy.ones(3) + 1.0))


def test_sa_h5py_13():
    """
    Test gridding tracks.
    """
    tracks = {"x" : numpy.array([10,20,30]),
              "y" : numpy.array([10,10,10]),
              "z" : numpy.array([-0.2,0.0,0.2])}

    h5_name = storm_analysis.getPathOutputTest("test_sa_hdf5.hdf5")
    
    # Tracks.
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(40, 40, 1, "")
        h5.addTracks(tracks)

    with saH5Py.SAH5Grid(filename = h5_name, scale = 1, z_bins = 3) as h5g:
        im_2d = h5g.gridTracks2D()
        im_3d = h5g.gridTracks3D(-0.201,0.201)

        for i in range(tracks["x"].size):
            assert(im_2d[int(tracks["x"][i]), int(tracks["y"][i])] == 1)
            assert(im_3d[int(tracks["x"][i]), int(tracks["y"][i]), i] == 1)

        im_3d = h5g.gridTracks3D(-0.1,0.1)
        assert(im_3d[int(tracks["x"][0]), int(tracks["y"][0]), 0] == 0)
        assert(im_3d[int(tracks["x"][1]), int(tracks["y"][1]), 1] == 1)
        assert(im_3d[int(tracks["x"][2]), int(tracks["y"][2]), 2] == 0)


def test_sa_h5py_14():
    """
    Test gridding tracks with dx, dy.
    """
    tracks = {"x" : numpy.array([10,20,30]),
              "y" : numpy.array([10,10,10]),
              "z" : numpy.array([-0.2,0.0,0.2])}

    dx = 1
    dy = 2
    h5_name = storm_analysis.getPathOutputTest("test_sa_hdf5.hdf5")
    
    # Tracks.
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(40, 40, 1, "")
        h5.addTracks(tracks)

    with saH5Py.SAH5Grid(filename = h5_name, scale = 1, z_bins = 3) as h5g:
        im_2d = h5g.gridTracks2D(dx = dx, dy = dy)
        im_3d = h5g.gridTracks3D(-0.201,0.201, dx = dx, dy = dy)

        for i in range(tracks["x"].size):
            assert(im_2d[int(tracks["x"][i]+dx), int(tracks["y"][i]+dy)] == 1)
            assert(im_3d[int(tracks["x"][i]+dx), int(tracks["y"][i]+dy), i] == 1)
        

def test_sa_h5py_15():
    """
    Test isAnalyzed()
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    empty = {"x" : numpy.array([]),
             "y" : numpy.array([])}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(100, 100, 2, "")
        h5.addLocalizations(peaks, 0)
        h5.addLocalizations(empty, 1)

    # Read data.
    with saH5Py.SAH5Py(h5_name) as h5:

        for i, elt in enumerate([True, True, False]):
            assert (h5.isAnalyzed(i) == elt)


def test_sa_h5py_16():
    """
    Test that localizations iterator skips empty frames.
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    empty = {"x" : numpy.array([]),
             "y" : numpy.array([])}    

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(100, 100, 5, "")
        h5.addLocalizations(peaks, 0)
        h5.addLocalizations(empty, 1)
        h5.addLocalizations(peaks, 2)

    # Read data.
    with saH5Py.SAH5Py(h5_name) as h5:
        for fnum, locs in h5.localizationsIterator():
            assert(fnum != 1)


def test_sa_h5py_17():
    """
    Test that localizations iterator does not skip empty frames (when requested not to).
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    empty = {"x" : numpy.array([]),
             "y" : numpy.array([])}    

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(100, 100, 5, "")
        h5.addLocalizations(peaks, 0)
        h5.addLocalizations(empty, 1)
        h5.addLocalizations(peaks, 2)

    # Read data.
    with saH5Py.SAH5Py(h5_name) as h5:
        for i, [fnum, locs] in enumerate(h5.localizationsIterator(skip_empty = False)):
            assert(i == fnum)


def test_sa_h5py_18():
    """
    Analysis finished flag.
    """
    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(100, 100, 5, "")
        assert(not h5.isAnalysisFinished())
        h5.setAnalysisFinished(True)
        assert(h5.isAnalysisFinished())
        h5.setAnalysisFinished(False)
        assert(not h5.isAnalysisFinished())


def test_sa_h5py_19():
    """
    Test getting specific fields.
    """
    peaks = {"bar" : numpy.zeros(10),
             "x" : numpy.zeros(10),
             "y" : numpy.zeros(10)}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)
    
    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(100, 100, 1, "")
        h5.addLocalizations(peaks, 0)

    # Get data.
    with saH5Py.SAH5Py(h5_name) as h5:
        locs = h5.getLocalizationsInFrame(0)
        for elt in ["bar", "x", "y"]:
            assert elt in locs

        locs = h5.getLocalizationsInFrame(0, fields = ["x"])
        assert "x" in locs
        for elt in ["bar", "y"]:
            assert not elt in locs
    
        
if (__name__ == "__main__"):
    test_sa_h5py_1()
    test_sa_h5py_2()
    test_sa_h5py_3()
    test_sa_h5py_4()
    test_sa_h5py_5()
    test_sa_h5py_6()
    test_sa_h5py_7()
    test_sa_h5py_8()
    test_sa_h5py_9()
    test_sa_h5py_10()
    test_sa_h5py_11()
    test_sa_h5py_12()
    test_sa_h5py_13()
    test_sa_h5py_14()
    test_sa_h5py_15()
    test_sa_h5py_16()
    test_sa_h5py_17()
    test_sa_h5py_18()
    test_sa_h5py_19()
    
