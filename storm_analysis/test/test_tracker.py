#!/usr/bin/env python
"""
Tests of sa_utilities.tracker
"""
import numpy

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_utilities.tracker as tracker


class FakeReader(object):
    """
    Fake movie reader object.
    """
    def __init__(self, n_frames = None, **kwds):
        super(FakeReader, self).__init__(**kwds)
        self.n_frames = n_frames

    def hashID(self):
        return "xyzzy"

    def getMovieL(self):
        return 10

    def getMovieX(self):
        return 128

    def getMovieY(self):
        return 128
    

def test_tracker_1():
    """
    Basic tracking test.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0]),
             "sum" : numpy.array([4.0, 4.0, 4.0])}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            temp = {}
            for elt in peaks:
                temp[elt] = peaks[elt][i:]
            h5.addLocalizations(temp, i)

        h5.addMovieInformation(FakeReader(n_frames = 3))

    # Track.
    tracker.tracker(h5_name, radius = 0.1)

    # Tracking.
    with saH5Py.SAH5Py(h5_name) as h5:

        # Check that we have the right number of tracks.
        assert(h5.getNTracks() == 3)

        # Check tracks.
        for t in h5.tracksIterator():
            assert(numpy.allclose(peaks["x"], t["x"]))
            assert(numpy.allclose(peaks["y"], t["y"]))
            assert(numpy.allclose(numpy.array([0,0,0]), t["frame_number"]))
            assert(numpy.allclose(numpy.array([1,2,3]), t["track_length"]))

        # Check that the localizations 'track_id' field is correct.
        for fnum, locs in h5.localizationsIterator(fields = ["track_id"]):
            assert(numpy.allclose(numpy.array([0,1,2])[fnum:], locs["track_id"]))


def test_tracker_2():
    """
    Test descriptor.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0]),
             "sum" : numpy.array([4.0, 4.0, 4.0])}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(4):
            temp = {}
            for elt in peaks:
                temp[elt] = peaks[elt][:i]
            if(len(temp["x"])>0):
                h5.addLocalizations(temp, i)

        h5.addMovieInformation(FakeReader(n_frames = 4))

    # Track.
    tracker.tracker(h5_name, descriptor = "1212", radius = 0.1)

    # Tracking.
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.getNTracks() == 3)
        for t in h5.tracksIterator():
            assert(numpy.allclose(peaks["x"], t["x"]))
            assert(numpy.allclose(peaks["y"], t["y"]))
            assert(numpy.allclose(numpy.array([1,0,1]), t["category"]))

        # Check localization categories.
        for fnum, locs in h5.localizationsIterator(fields = ["category"]):
            if (fnum == 1):
                assert(numpy.allclose(numpy.array([1]), locs["category"]))
            if (fnum == 2):
                assert(numpy.allclose(numpy.array([0, 0]), locs["category"]))
            if (fnum == 3):
                assert(numpy.allclose(numpy.array([1, 1, 1]), locs["category"]))


def test_tracker_3():
    """
    Test tracking over an activation frame.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0]),
             "sum" : numpy.array([4.0, 4.0, 4.0])}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addLocalizations(peaks, 0)
        h5.addLocalizations(peaks, 2)
        h5.addMovieInformation(FakeReader(n_frames = 3))

    # Track.
    tracker.tracker(h5_name, descriptor = "101", radius = 0.1)

    # Tracking.
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.getNTracks() == 3)
        for t in h5.tracksIterator():
            assert(numpy.allclose(peaks["x"], t["x"]))
            assert(numpy.allclose(peaks["y"], t["y"]))
            assert(numpy.allclose(numpy.array([2,2,2]), t["track_length"]))
            

def test_tracker_4():
    """
    Test that nearest object is assigned to the nearest track.
    """
    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
                 "y" : numpy.array([1.0, 1.0, 1.0]),
                 "sum" : numpy.array([4.0, 4.0, 4.0])}
        h5.addLocalizations(peaks, 0)
        
        peaks = {"x" : numpy.array([1.0, 2.0, 2.1, 3.0]),
                 "y" : numpy.array([1.0, 1.0, 1.0, 1.0]),
                 "sum" : numpy.array([4.0, 4.0, 5.0, 4.0])}
        h5.addLocalizations(peaks, 1)

        h5.addMovieInformation(FakeReader(n_frames = 2))

    # Track.
    tracker.tracker(h5_name, radius = 0.5)

    # Tracking.
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.getNTracks() == 4)
        for t in h5.tracksIterator():
            assert(numpy.allclose(numpy.array([1,2,3,2.1]), t["x"]))
            assert(numpy.allclose(numpy.array([2,2,2,1]), t["track_length"]))


def test_tracker_5():
    """
    Test that nearest track is assigned to the nearest object.
    """
    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
                 "y" : numpy.array([1.0, 1.0, 1.0]),
                 "sum" : numpy.array([4.0, 4.0, 4.0])}
        h5.addLocalizations(peaks, 0)
        
        peaks = {"x" : numpy.array([2.0]),
                 "y" : numpy.array([1.0]),
                 "sum" : numpy.array([4.0])}
        h5.addLocalizations(peaks, 1)

        h5.addMovieInformation(FakeReader(n_frames = 2))

    # Track.
    tracker.tracker(h5_name, radius = 1.1)

    # Tracking.
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.getNTracks() == 3)
        for t in h5.tracksIterator():
            assert(numpy.allclose(numpy.array([1,3,2]), t["x"]))
            assert(numpy.allclose(numpy.array([0,2,1]), t["track_id"]))
            assert(numpy.allclose(numpy.array([1,1,2]), t["track_length"]))


def test_tracker_6():
    """
    Test max_gap parameter.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0]),
             "sum" : numpy.array([4.0, 4.0, 4.0])}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addLocalizations(peaks, 0)
        h5.addLocalizations(peaks, 2)
        h5.addMovieInformation(FakeReader(n_frames = 3))

    # Track.
    tracker.tracker(h5_name, radius = 0.1)

    # Tracking.
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.getNTracks() == 6)
        for t in h5.tracksIterator():
            assert(numpy.allclose(t["track_length"], numpy.ones(6)))

    # Redo the tracking allowing single frame gaps.
    tracker.tracker(h5_name, max_gap = 1, radius = 0.1)

    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.getNTracks() == 3)
        for t in h5.tracksIterator():
            assert(numpy.allclose(t["track_length"], 2.0*numpy.ones(3)))


def test_tracker_7():
    """
    Test handling of 0.0 radius.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0]),
             "sum" : numpy.array([4.0, 4.0, 4.0])}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(4):
            temp = {}
            for elt in peaks:
                temp[elt] = peaks[elt][:i]
            if(len(temp["x"])>0):
                h5.addLocalizations(temp, i)

        h5.addMovieInformation(FakeReader(n_frames = 4))

    # Track.
    tracker.tracker(h5_name, descriptor = "1212", radius = 0.0)

    # Tracking.
    with saH5Py.SAH5Py(h5_name) as h5:

        # Check that there are no tracks.
        assert(h5.getNTracks() == 0)
        
        # Check localization categories.
        for fnum, locs in h5.localizationsIterator(fields = ["category"]):
            if (fnum == 1):
                assert(numpy.allclose(numpy.array([1]), locs["category"]))
            if (fnum == 2):
                assert(numpy.allclose(numpy.array([0, 0]), locs["category"]))
            if (fnum == 3):
                assert(numpy.allclose(numpy.array([1, 1, 1]), locs["category"]))


def test_tracker_8():
    """
    Test tracking over an empty frame.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0]),
             "sum" : numpy.array([4.0, 4.0, 4.0])}

    empty = {"x" : numpy.array([]),
             "y" : numpy.array([]),
             "sum" : numpy.array([])} 

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.addLocalizations(peaks, 0)
        h5.addLocalizations(empty, 1)
        h5.addLocalizations(peaks, 2)
        h5.addMovieInformation(FakeReader(n_frames = 3))

    # Track.
    tracker.tracker(h5_name, descriptor = "111", radius = 0.1)

    # Tracking.
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.getNTracks() == 6)
        for t in h5.tracksIterator():
            assert(numpy.allclose(numpy.ones(6), t["track_length"]))
            
            
if (__name__ == "__main__"):
    test_tracker_1()
    test_tracker_2()
    test_tracker_3()
    test_tracker_4()
    test_tracker_5()
    test_tracker_6()
    test_tracker_7()
    test_tracker_8()

    
