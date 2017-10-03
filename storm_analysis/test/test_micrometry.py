#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.writeinsight3 as writeinsight3

import storm_analysis.micrometry.micrometry as micrometry


def test_micrometry_1():
    """
    Test micrometry on matching data.
    """

    locs1_name = storm_analysis.getPathOutputTest("locs1.bin")
    locs2_name = storm_analysis.getPathOutputTest("locs2.bin")

    # Create test data.
    im_size = 512
    n_points = 50

    numpy.random.seed(0)

    m_data = i3dtype.createDefaultI3Data(n_points)

    i3dtype.posSet(m_data, "x", numpy.random.uniform(high = im_size, size = n_points))
    i3dtype.posSet(m_data, "y", numpy.random.uniform(high = im_size, size = n_points))

    with writeinsight3.I3Writer(locs1_name) as i3w:
        i3w.addMolecules(m_data)
    with writeinsight3.I3Writer(locs2_name) as i3w:
        i3w.addMolecules(m_data)

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

    locs1_name = storm_analysis.getPathOutputTest("locs1.bin")
    locs2_name = storm_analysis.getPathOutputTest("locs2.bin")

    # Create test data.
    im_size = 512
    n_points = 50

    numpy.random.seed(0)

    with writeinsight3.I3Writer(locs1_name) as i3w:
        m_data = i3dtype.createDefaultI3Data(n_points)
        i3dtype.posSet(m_data, "x", numpy.random.uniform(high = im_size, size = n_points))
        i3dtype.posSet(m_data, "y", numpy.random.uniform(high = im_size, size = n_points))
        i3w.addMolecules(m_data)
        
    with writeinsight3.I3Writer(locs2_name) as i3w:
        m_data = i3dtype.createDefaultI3Data(n_points)
        i3dtype.posSet(m_data, "x", numpy.random.uniform(high = im_size, size = n_points))
        i3dtype.posSet(m_data, "y", numpy.random.uniform(high = im_size, size = n_points))
        i3w.addMolecules(m_data)

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
