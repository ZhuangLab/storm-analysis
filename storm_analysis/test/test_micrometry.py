#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.micrometry.micrometry as micrometry


def test_micrometry_1():
    """
    Test micrometry on matching data.
    """

    locs1_name = storm_analysis.getPathOutputTest("locs1.hdf5")
    locs2_name = storm_analysis.getPathOutputTest("locs2.hdf5")

    # Create test data.
    im_size = 512
    n_points = 50

    numpy.random.seed(0)

    locs = {"x" : numpy.random.uniform(high = im_size, size = n_points),
            "y" : numpy.random.uniform(high = im_size, size = n_points)}
    
    with saH5Py.SAH5Py(locs1_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(512, 512, 1, "")
        h5.addLocalizations(locs, 0)
    with saH5Py.SAH5Py(locs2_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(512, 512, 1, "")
        h5.addLocalizations(locs, 0)

    # Test
    mm = micrometry.Micrometry(locs1_name,
                               min_size = 5.0,
                               max_size = 100.0,
                               max_neighbors = 20)
    [best_ratio, best_transform] = mm.findTransform(locs2_name, 1.0e-2)

    assert(best_ratio > 10.0)
    expected = [[0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0]]
    for i, elt in enumerate(best_transform):
        for j in range(3):
            assert(abs(expected[i][j] - best_transform[i][j]) < 1.0e-6)
    

def test_micrometry_2():
    """
    Test micrometry on random data.
    """

    locs1_name = storm_analysis.getPathOutputTest("locs1.hdf5")
    locs2_name = storm_analysis.getPathOutputTest("locs2.hdf5")

    # Create test data.
    im_size = 512
    n_points = 50

    numpy.random.seed(0)

    with saH5Py.SAH5Py(locs1_name, is_existing = False, overwrite = True) as h5:
        locs = {"x" : numpy.random.uniform(high = im_size, size = n_points),
                "y" : numpy.random.uniform(high = im_size, size = n_points)}
        h5.setMovieInformation(512, 512, 1, "")
        h5.addLocalizations(locs, 0)

    with saH5Py.SAH5Py(locs2_name, is_existing = False, overwrite = True) as h5:
        locs = {"x" : numpy.random.uniform(high = im_size, size = n_points),
                "y" : numpy.random.uniform(high = im_size, size = n_points)}
        h5.setMovieInformation(512, 512, 1, "")
        h5.addLocalizations(locs, 0)

    # Test
    mm = micrometry.Micrometry(locs1_name,
                               min_size = 5.0,
                               max_size = 100.0,
                               max_neighbors = 20)
    [best_ratio, best_transform] = mm.findTransform(locs2_name, 1.0e-2)

    assert(best_ratio < 10.0)
            
            
if (__name__ == "__main__"):
    test_micrometry_1()
    test_micrometry_2()
