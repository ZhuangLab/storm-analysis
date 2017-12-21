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
    Simple tracking test 1.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0]),
             "sum" : numpy.array([4.0, 4.0, 4.0])}

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name) as h5:
        for i in range(3):
            temp = {}
            for elt in peaks:
                temp[elt] = peaks[elt][i:]
            print(i, temp["x"])
            h5.addLocalizations(temp, i)

        h5.addMovieInformation(FakeReader(n_frames = 3))

    # Track.
    tracker.tracker(h5_name, radius = 0.1)

    # Tracking.
    with saH5Py.SAH5Py(h5_name) as h5:
        print(h5.getNTracks())
        #assert(h5.getNTracks() == 3)
        for tracks in h5.tracksIterator():
            print(tracks)
        

if (__name__ == "__main__"):
    test_tracker_1()
    
