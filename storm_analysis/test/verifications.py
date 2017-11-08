#!/usr/bin/env python
"""
These are used at during testing to check the results.

Hazen 02/17
"""
import numpy

import storm_analysis.sa_library.readinsight3 as readinsight3


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

def verifyNumberLocalizations(bin_fname):
    """
    Return the number of localizations in a I3 file.
    """
    return readinsight3.loadI3File(bin_fname, verbose = False).size

def verifyZWasCalculated(bin_fname):
    """
    Return the true if all the Z values are not exactly identical.
    """
    z_vals = readinsight3.loadI3File(bin_fname, verbose = False)["z"]
    return (numpy.std(z_vals) > 1.0e-6)


