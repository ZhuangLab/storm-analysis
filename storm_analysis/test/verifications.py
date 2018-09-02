#!/usr/bin/env python
"""
These are used at during testing to check the results.

Hazen 02/17
"""
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py


def verifyDriftCorrection(actual_drift_fname, measured_drift_fname):
    """
    Return maximum difference between drift correction files.
    """

    actual_drift = numpy.loadtxt(actual_drift_fname)
    measured_drift = numpy.loadtxt(measured_drift_fname)

    diffs = []
    for i in range(4):
        diffs.append(numpy.max(numpy.abs(actual_drift[:,i] - measured_drift[:,i])))

    return diffs

def verifyIsCloseEnough(number1, number2, margin = 0.05):
    """
    Return true if number1 is within margin of number 2.
    """
    max_diff = number2 * margin
    return (abs(number1 - number2) < max_diff)

def verifyNumberLocalizations(h5_name):
    """
    Return the number of localizations in a HDF5 file.
    """
    n_locs = None
    with saH5Py.SAH5Py(h5_name) as h5:
        assert(h5.isAnalysisFinished())
        n_locs = h5.getNLocalizations()
    return n_locs

def verifyZWasCalculated(h5_name):
    """
    Return the true if all the Z values are not exactly identical.
    """
    locs = None
    with saH5Py.SAH5Py(h5_name) as h5:
        locs = h5.getLocalizations(fields = ["z"])
    return (numpy.std(locs["z"]) > 1.0e-6)


