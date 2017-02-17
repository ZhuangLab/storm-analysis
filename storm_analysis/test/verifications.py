#!/usr/bin/env python
"""
Checks drift correction results for accuracy.

Hazen 02/17
"""
import numpy

import storm_analysis.sa_library.readinsight3 as readinsight3

def verifyDriftCorrection(actual_drift_fname, measured_drift_fname):

    actual_drift = numpy.loadtxt(actual_drift_fname)
    measured_drift = numpy.loadtxt(measured_drift_fname)

    diffs = []
    for i in range(4):
        diffs.append(numpy.max(numpy.abs(actual_drift[:,i] - measured_drift[:,i])))

    return diffs

def verifyNumberLocalizations(bin_fname):
    return readinsight3.loadI3File(bin_fname, verbose = False).size
