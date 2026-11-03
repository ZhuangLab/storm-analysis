#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.micrometry.micrometry as micrometry
import storm_analysis.micrometry.quads as quads


def similarityTransform(rotation, scale, shift, mirror):
    """
    Returns the 2x2 matrix and the shift of a similarity transform, with
    rotation in degrees, mirroring x first if mirror is True.
    """
    angle = numpy.deg2rad(rotation)
    matrix = scale * numpy.array([[numpy.cos(angle), -numpy.sin(angle)],
                                  [numpy.sin(angle), numpy.cos(angle)]])
    if mirror:
        matrix = numpy.matmul(matrix, numpy.diag([-1.0, 1.0]))
    return [matrix, numpy.array(shift)]


def transformError(transform, matrix, shift, im_size = 512):
    """
    The RMS distance across the field between where the 'reference' to
    'other' part of transform puts a point and where it should be.
    """
    grid = numpy.stack(numpy.meshgrid(numpy.linspace(0, im_size, 11), numpy.linspace(0, im_size, 11)), axis = -1).reshape(-1, 2)
    estimate = micrometry.transformPoints(grid, transform[2:])
    expected = numpy.matmul(grid, matrix.T) + shift
    return numpy.sqrt(numpy.mean(numpy.sum((estimate - expected)**2, axis = 1)))


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


def test_micrometry_3():
    """
    Test that makeQuads() gives the same codes as makeQuad(), and the
    same groups of points for a mirrored copy of the points.
    """
    rng = numpy.random.default_rng(1)
    points = rng.uniform(high = 512, size = (100, 2))

    [indices, codes] = quads.makeQuads(points)
    assert(indices.shape[0] > 1000)
    for i in range(indices.shape[0]):
        code = quads.makeQuad(*points[indices[i]])
        assert code is not None
        assert numpy.allclose(code, codes[i], atol = 1.0e-12)

    # Each group of four points appears once.
    groups = set(tuple(sorted(elt)) for elt in indices.tolist())
    assert(len(groups) == indices.shape[0])

    # A mirror image has the same groups, with different codes.
    [m_indices, m_codes] = quads.makeQuads(points * numpy.array([-1.0, 1.0]))
    assert(groups == set(tuple(sorted(elt)) for elt in m_indices.tolist()))


def test_micrometry_4():
    """
    Test micrometry on mirrored, rotated and scaled data with
    localization error and missing localizations.
    """
    rng = numpy.random.default_rng(2)
    im_size = 512
    n_points = 100
    sigma = 0.5

    [matrix, shift] = similarityTransform(8.0, 1.2, [400.0, -50.0], True)
    points = rng.uniform(high = im_size, size = (n_points, 2))
    other = numpy.matmul(points, matrix.T) + shift

    # 30% of the beads are missing from each list.
    keep1 = (rng.uniform(size = n_points) > 0.3)
    keep2 = (rng.uniform(size = n_points) > 0.3)
    points1 = points[keep1] + rng.normal(scale = sigma, size = (numpy.count_nonzero(keep1), 2))
    points2 = other[keep2] + rng.normal(scale = sigma, size = (numpy.count_nonzero(keep2), 2))
    points2 = rng.permutation(points2)

    [score, transform] = micrometry.findTransformFromPoints(points1, points2, localization_error = sigma)

    assert(score > 10.0)
    assert(transformError(transform, matrix, shift) < 0.5)

    # The 'other' to 'reference' transform is the inverse.
    round_trip = micrometry.transformPoints(micrometry.transformPoints(points1, transform[2:]), transform[:2])
    assert(numpy.max(numpy.abs(round_trip - points1)) < 0.1)


def test_micrometry_5():
    """
    Test that unrelated lists are rejected for a range of field sizes,
    localization counts and localization errors.
    """
    for [n_points, im_size, localization_error, seed] in [[20, 512, 0.5, 3],
                                                          [50, 256, 0.5, 4],
                                                          [50, 512, 1.0, 5],
                                                          [200, 512, 0.1, 6],
                                                          [400, 2048, 0.5, 7]]:
        rng = numpy.random.default_rng(seed)
        points1 = rng.uniform(high = im_size, size = (n_points, 2))
        points2 = rng.uniform(high = im_size, size = (n_points, 2))
        [score, transform] = micrometry.findTransformFromPoints(points1, points2, localization_error = localization_error)
        assert(score < 10.0)

    
if (__name__ == "__main__"):
    test_micrometry_1()
    test_micrometry_2()
    test_micrometry_3()
    test_micrometry_4()
    test_micrometry_5()
