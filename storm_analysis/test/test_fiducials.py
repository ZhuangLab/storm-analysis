#!/usr/bin/env python
"""
Tests of sa_utilities.fiducials
"""
import numpy

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_utilities.fiducials as fiducials


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
        return self.n_frames

    def getMovieX(self):
        return 128

    def getMovieY(self):
        return 128
    

def test_fiducials_1():
    """
    Basic fiducials test.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            h5.addLocalizations(peaks, i)

        h5.addMovieInformation(FakeReader(n_frames = 3))
        
    # Track fiducials..
    fiducials.trackFiducials(h5_name, radius = 0.1)

    # Check.
    with saH5Py.SAH5Py(h5_name) as h5:
        for fnum, locs in h5.localizationsIterator(fields = ["fiducial_id"]):
            assert(numpy.allclose(locs["fiducial_id"], numpy.arange(3)))


def test_fiducials_2():
    """
    Basic fiducials test.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            temp = {}
            for elt in peaks:
                temp[elt] = peaks[elt][i:]
            h5.addLocalizations(temp, i)

        h5.addMovieInformation(FakeReader(n_frames = 4))
        
    # Track fiducials..
    fiducials.trackFiducials(h5_name, radius = 0.1)

    # Check.
    with saH5Py.SAH5Py(h5_name) as h5:
        expected = numpy.array([0,1,2])
        for fnum, locs in h5.localizationsIterator(fields = ["fiducial_id"]):
            assert numpy.allclose(locs["fiducial_id"], expected[fnum:])
            

def test_fiducials_3():
    """
    Basic fiducials test.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            temp = {}
            for elt in peaks:
                temp[elt] = peaks[elt][i:]
            h5.addLocalizations(temp, i)

        h5.addMovieInformation(FakeReader(n_frames = 4))
        
    # Track fiducials..
    fiducials.trackFiducials(h5_name, radius = 0.1, reference_frame = 2)

    # Check.
    with saH5Py.SAH5Py(h5_name) as h5:
        expected = numpy.array([-1,-1,0])
        for fnum, locs in h5.localizationsIterator(fields = ["fiducial_id"]):
            assert numpy.allclose(locs["fiducial_id"], expected[fnum:])


def test_fiducials_4():
    """
    Test no localizations in reference frame.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            h5.addLocalizations(peaks, i)
        h5.addMovieInformation(FakeReader(n_frames = 5))
        
    # Track fiducials..
    okay = False
    try:
        fiducials.trackFiducials(h5_name, radius = 0.1, reference_frame = 3)
    except fiducials.FiducialException:
        okay = True
    assert okay


def test_fiducials_5():
    """
    Iterator test.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            h5.addLocalizations(peaks, i)

        h5.addMovieInformation(FakeReader(n_frames = 6))
        
    # Track fiducials..
    fiducials.trackFiducials(h5_name, radius = 0.1)

    # Check.
    with fiducials.SAH5Fiducials(h5_name) as h5:
        for fdcl in h5.fiducialsIterator():
            assert(numpy.allclose(fdcl["frame"], numpy.arange(3)))


def test_fiducials_6():
    """
    Iterator test.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            temp = {}
            for elt in peaks:
                temp[elt] = peaks[elt][i:]
            h5.addLocalizations(temp, i)

        h5.addMovieInformation(FakeReader(n_frames = 4))
        
    # Track fiducials..
    fiducials.trackFiducials(h5_name, radius = 0.1)

    # Check.
    with fiducials.SAH5Fiducials(h5_name) as h5:
        expected = numpy.array([0,1,2])
        for i, fdcl in enumerate(h5.fiducialsIterator()):
            assert(numpy.allclose(fdcl["frame"], numpy.arange(i+1)))

            
def test_fiducials_7():
    """
    Iterator test.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            temp = {}
            for elt in peaks:
                temp[elt] = peaks[elt][i:]
            h5.addLocalizations(temp, i)

        h5.addMovieInformation(FakeReader(n_frames = 4))
        
    # Track fiducials..
    fiducials.trackFiducials(h5_name, radius = 0.1, reference_frame = 2)

    # Check.
    with fiducials.SAH5Fiducials(h5_name) as h5:
        for fdcl in h5.fiducialsIterator():
            assert(numpy.allclose(fdcl["frame"], numpy.arange(3)))


def test_fiducials_8():
    """
    Gap test.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in [0,1,3]:
            h5.addLocalizations(peaks, i)

        h5.addMovieInformation(FakeReader(n_frames = 4))
        
    # Track fiducials..
    fiducials.trackFiducials(h5_name, radius = 0.1, max_gap = 1)

    # Check.
    with fiducials.SAH5Fiducials(h5_name) as h5:
        expected = numpy.array([0,1,3])
        for fdcl in h5.fiducialsIterator():
            assert(numpy.allclose(fdcl["frame"], expected))
            

def test_fiducials_9():
    """
    Test fiducial averaging.
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            h5.addLocalizations(peaks, i)

        h5.addMovieInformation(FakeReader(n_frames = 3))
        
    # Track fiducials..
    fiducials.trackFiducials(h5_name, radius = 0.1)

    # Check
    with fiducials.SAH5Fiducials(h5_name) as h5:
        [ave, n] = h5.averageFiducials(fields = ["y"])
        assert(numpy.allclose(ave["y"], numpy.ones(3)))


def test_fiducials_10():
    """
    Test fiducial averaging (preload_all = False).
    """
    peaks = {"x" : numpy.array([1.0, 2.0, 3.0]),
             "y" : numpy.array([1.0, 1.0, 1.0])}

    filename = "test_fiducials.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write data.
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        for i in range(3):
            h5.addLocalizations(peaks, i)

        h5.addMovieInformation(FakeReader(n_frames = 3))
        
    # Track fiducials..
    fiducials.trackFiducials(h5_name, radius = 0.1)

    # Check
    with fiducials.SAH5Fiducials(h5_name) as h5:
        [ave, n] = h5.averageFiducials(fields = ["y"], preload_all = False)
        assert(numpy.allclose(ave["y"], numpy.ones(3)))
    
if (__name__ == "__main__"):
    test_fiducials_1()
    test_fiducials_2()
    test_fiducials_3()
    test_fiducials_4()
    test_fiducials_5()
    test_fiducials_6()
    test_fiducials_7()
    test_fiducials_8()
    test_fiducials_9()
    test_fiducials_10()
