#!/usr/bin/env python
import numpy
import os
import pickle

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.micrometry.make_test_data as makeTestData
import storm_analysis.micrometry.micrometry as micrometry
import storm_analysis.micrometry.quads as quads
import storm_analysis.micrometry.refine_map as refineMap


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


def test_micrometry_6():
    """
    Test runMicrometry() on HDF5 files, and that it saves the mapping.
    """
    locs1_name = storm_analysis.getPathOutputTest("locs1.hdf5")
    locs2_name = storm_analysis.getPathOutputTest("locs2.hdf5")
    map_name = storm_analysis.getPathOutputTest("micrometry.map")

    im_size = 512
    makeTestData.makeTestData(locs1_name, locs2_name, 80, 0.8,
                              im_size = im_size,
                              rotation = -5.0,
                              scale = 0.9,
                              dx = 10.0,
                              dy = 5.0,
                              mirror = True,
                              localization_error = 0.2,
                              seed = 8)

    [score, transform] = micrometry.runMicrometry(locs1_name, locs2_name, map_name, no_plots = True)
    assert(score > 10.0)

    # makeTestData() transforms about the center of the field.
    [matrix, shift] = similarityTransform(-5.0, 0.9, [10.0, 5.0], True)
    center = numpy.array([0.5 * im_size, 0.5 * im_size])
    shift = shift + center - numpy.matmul(matrix, center)
    assert(transformError(transform, matrix, shift) < 0.5)

    with open(map_name, 'rb') as fp:
        mapping = pickle.load(fp)
    for [i, key] in enumerate(["1_0_x", "1_0_y", "0_1_x", "0_1_y"]):
        assert numpy.allclose(mapping[key], transform[i])

    # The default Micrometry() and the array entry point give the same answer.
    mm = micrometry.Micrometry(locs1_name, verbose = False)
    [h5_score, h5_transform] = mm.findTransform(locs2_name)
    [array_score, array_transform] = micrometry.findTransformFromPoints(micrometry.loadLocalizations(locs1_name),
                                                                        micrometry.loadLocalizations(locs2_name))
    assert(abs(h5_score - array_score) < 1.0e-6)
    assert numpy.allclose(h5_transform, array_transform)


def test_micrometry_7():
    """
    Test refining a mapping that is a few pixels off.
    """
    locs1_name = storm_analysis.getPathOutputTest("locs1.hdf5")
    locs2_name = storm_analysis.getPathOutputTest("locs2.hdf5")
    map_name = storm_analysis.getPathOutputTest("micrometry.map")
    refined_name = storm_analysis.getPathOutputTest("micrometry_refined.map")

    makeTestData.makeTestData(locs1_name, locs2_name, 100, 0.9,
                              rotation = 3.0,
                              dx = -20.0,
                              localization_error = 0.2,
                              seed = 9)

    # A mapping that is right but for a 1.5 pixel shift.
    [matrix, shift] = similarityTransform(3.0, 1.0, [-20.0, 0.0], False)
    center = numpy.array([256.0, 256.0])
    shift = shift + center - numpy.matmul(matrix, center)
    inverse = numpy.linalg.inv(matrix)
    inverse_shift = -numpy.matmul(inverse, shift) + numpy.array([1.5, 0.0])
    micrometry.saveMapping([numpy.array([inverse_shift[0], inverse[0,0], inverse[0,1]]),
                            numpy.array([inverse_shift[1], inverse[1,0], inverse[1,1]]),
                            numpy.array([shift[0], matrix[0,0], matrix[0,1]]),
                            numpy.array([shift[1], matrix[1,0], matrix[1,1]])],
                           map_name)

    n_pairs = refineMap.refineMap(locs1_name, locs2_name, map_name, refined_name)
    assert(n_pairs > 70)

    with open(refined_name, 'rb') as fp:
        mapping = pickle.load(fp)
    transform = [mapping["1_0_x"], mapping["1_0_y"], mapping["0_1_x"], mapping["0_1_y"]]
    assert(transformError(transform, matrix, shift) < 0.2)


def test_micrometry_8():
    """
    Test the estimated mapping error, which runMicrometry() warns about
    when it is more than a third of the pairing radius.
    """
    rng = numpy.random.default_rng(10)
    sigma = 0.5
    limit = micrometry.pairingRadius(sigma)/3.0
    [matrix, shift] = similarityTransform(5.0, 1.1, [20.0, -10.0], False)

    # Beads spread over the field, all seen in both channels.
    points = rng.uniform(high = 512, size = (60, 2))
    mm = micrometry.Micrometry(ref_points = points + rng.normal(scale = sigma, size = points.shape), verbose = False)
    [score, transform] = mm.matchPoints(numpy.matmul(points, matrix.T) + shift + rng.normal(scale = sigma, size = points.shape),
                                        localization_error = sigma)
    assert(score > 10.0)
    assert(mm.getPairCount() == 60)
    assert(mm.getMappingError() < limit)

    # The second channel only has beads in one corner of the field.
    corner = rng.uniform(high = 100, size = (12, 2))
    points1 = numpy.concatenate((rng.uniform(high = 512, size = (40, 2)), corner))
    points1 = points1 + rng.normal(scale = sigma, size = points1.shape)
    points2 = numpy.matmul(corner, matrix.T) + shift + rng.normal(scale = sigma, size = corner.shape)

    mm = micrometry.Micrometry(ref_points = points1, verbose = False)
    [score, transform] = mm.matchPoints(points2, localization_error = sigma)
    assert(score > 10.0)
    assert(mm.getPairCount() == 12)
    assert(mm.getMappingError() > limit)

    # runMicrometry() warns, and still saves the mapping.
    locs1_name = storm_analysis.getPathOutputTest("locs1.hdf5")
    locs2_name = storm_analysis.getPathOutputTest("locs2.hdf5")
    map_name = storm_analysis.getPathOutputTest("micrometry_corner.map")
    for [filename, points] in [[locs1_name, points1], [locs2_name, points2]]:
        with saH5Py.SAH5Py(filename, is_existing = False, overwrite = True) as h5:
            h5.setMovieInformation(512, 512, 1, "")
            h5.addLocalizations({"x" : points[:,0], "y" : points[:,1]}, 0)
    if os.path.exists(map_name):
        os.remove(map_name)

    [score, transform] = micrometry.runMicrometry(locs1_name, locs2_name, map_name, localization_error = sigma, no_plots = True)
    assert(score > 10.0)
    assert os.path.exists(map_name)


if (__name__ == "__main__"):
    test_micrometry_1()
    test_micrometry_2()
    test_micrometry_3()
    test_micrometry_4()
    test_micrometry_5()
    test_micrometry_6()
    test_micrometry_7()
    test_micrometry_8()
