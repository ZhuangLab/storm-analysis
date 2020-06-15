#!/usr/bin/env python
"""
Test noise and recall fraction calculcations.
"""
import h5py
import numpy

import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_utilities.recall_fraction as rfrac


eps = 1.0e-6


def test_recall_1():
    """
    Test recall fraction calculation.
    """
    px = numpy.arange(10)
    py = numpy.ones(10)

    # Write GT data.
    gt_name = storm_analysis.getPathOutputTest("trf_gt")
    with saH5Py.SAH5Py(gt_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(1,1,1,"")
        h5.addLocalizations({"x" : px, "y" : py}, 0)

    for i in range(2,10):

        # Write test data.
        meas_name = storm_analysis.getPathOutputTest("trf_meas")
        with saH5Py.SAH5Py(meas_name, is_existing = False, overwrite = True) as h5:
            h5.setMovieInformation(1,1,1,"")
            h5.addLocalizations({"x" : px[:i], "y" : py[:i] + eps}, 0)

        # Test recall calculation.
        with saH5Py.SAH5Reader(gt_name) as h5_gt:
            with saH5Py.SAH5Reader(meas_name) as h5_meas:
                [rl, tl] = rfrac.recallFraction(h5_gt, h5_meas, 0.1)
                assert(rl == i)
                assert(tl == 10)


def test_recall_2():
    """
    Test recall fraction calculation.
    """
    px = numpy.arange(10)
    py = numpy.ones(10)

    # Write GT data.
    gt_name = storm_analysis.getPathOutputTest("trf_gt")
    with saH5Py.SAH5Py(gt_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(1,1,1,"")
        h5.addLocalizations({"x" : px, "y" : py}, 0)

    for i in range(2,10):

        # Write test data.
        meas_name = storm_analysis.getPathOutputTest("trf_meas")
        with saH5Py.SAH5Py(meas_name, is_existing = False, overwrite = True) as h5:
            h5.setMovieInformation(1,1,1,"")
            ty = numpy.copy(py) + eps
            ty[i:] += 0.2
            h5.addLocalizations({"x" : px, "y" : ty}, 0)

        # Test recall calculation.
        with saH5Py.SAH5Reader(gt_name) as h5_gt:
            with saH5Py.SAH5Reader(meas_name) as h5_meas:
                [rl, tl] = rfrac.recallFraction(h5_gt, h5_meas, 0.1)
                assert(rl == i)
                assert(tl == 10)


def test_recall_3():
    """
    Test recall fraction calculation.
    """
    px = numpy.arange(10)
    py = numpy.ones(10)

    # Write GT data.
    gt_name = storm_analysis.getPathOutputTest("trf_gt")
    with saH5Py.SAH5Py(gt_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(1,1,1,"")
        h5.addLocalizations({"x" : px, "y" : py}, 0)

    # Write test data.
    meas_name = storm_analysis.getPathOutputTest("trf_meas")
    with saH5Py.SAH5Py(meas_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(1,1,1,"")

    # Test recall calculation.
    with saH5Py.SAH5Reader(gt_name) as h5_gt:
        with saH5Py.SAH5Reader(meas_name) as h5_meas:
            [rl, tl] = rfrac.recallFraction(h5_gt, h5_meas, 0.1)
            assert (rl == 0)
            assert (tl == 10)

                
def test_noise_1():
    """
    Test noise fraction calculation.
    """
    px = numpy.arange(10)
    py = numpy.ones(10)

    # Write GT data.
    gt_name = storm_analysis.getPathOutputTest("trf_gt")
    with saH5Py.SAH5Py(gt_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(1,1,1,"")
        h5.addLocalizations({"x" : px, "y" : py}, 0)

    for i in range(2,10):

        # Write test data.
        meas_name = storm_analysis.getPathOutputTest("trf_meas")
        with saH5Py.SAH5Py(meas_name, is_existing = False, overwrite = True) as h5:
            h5.setMovieInformation(1,1,1,"")
            h5.addLocalizations({"x" : px[:i], "y" : py[:i] + eps}, 0)

        # Test noise calculation.
        with saH5Py.SAH5Reader(gt_name) as h5_gt:
            with saH5Py.SAH5Reader(meas_name) as h5_meas:
                [nl, tl] = rfrac.noiseFraction(h5_gt, h5_meas, 0.1)
                assert(nl == 0)
                assert(tl == i)


def test_noise_2():
    """
    Test noise fraction calculation.
    """
    px = numpy.arange(10)
    py = numpy.ones(10)

    # Write GT data.
    gt_name = storm_analysis.getPathOutputTest("trf_gt")
    with saH5Py.SAH5Py(gt_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(1,1,1,"")
        h5.addLocalizations({"x" : px, "y" : py}, 0)

    for i in range(2,10):

        # Write test data.
        meas_name = storm_analysis.getPathOutputTest("trf_meas")
        with saH5Py.SAH5Py(meas_name, is_existing = False, overwrite = True) as h5:
            h5.setMovieInformation(1,1,1,"")
            ty = numpy.copy(py) + eps
            ty[:i] += 0.2
            h5.addLocalizations({"x" : px, "y" : ty}, 0)

        # Test noise calculation.
        with saH5Py.SAH5Reader(gt_name) as h5_gt:
            with saH5Py.SAH5Reader(meas_name) as h5_meas:
                [nl, tl] = rfrac.noiseFraction(h5_gt, h5_meas, 0.1)
                assert(nl == i)
                assert(tl == 10)


def test_noise_3():
    """
    Test noise fraction calculation.
    """
    px = numpy.arange(10)
    py = numpy.ones(10)

    # Write GT data.
    gt_name = storm_analysis.getPathOutputTest("trf_gt")
    with saH5Py.SAH5Py(gt_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(1,1,1,"")

    # Write test data.
    meas_name = storm_analysis.getPathOutputTest("trf_meas")
    with saH5Py.SAH5Py(meas_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(1,1,1,"")
        h5.addLocalizations({"x" : px, "y" : py}, 0)

    # Test noise calculation.
    with saH5Py.SAH5Reader(gt_name) as h5_gt:
        with saH5Py.SAH5Reader(meas_name) as h5_meas:
            [nl, tl] = rfrac.noiseFraction(h5_gt, h5_meas, 0.1)
            assert (nl == 10)
            assert (tl == 10)


if (__name__ == "__main__"):
    test_recall_1()
    test_recall_2()
    test_recall_3()
    test_noise_1()
    test_noise_2()
    test_noise_3()
