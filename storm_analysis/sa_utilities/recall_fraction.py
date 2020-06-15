#!/usr/bin/env python
"""
Functions for calculating recall, etc.

Hazen 09/17
"""
import numpy

import storm_analysis.sa_library.ia_utilities_c as iaUtilsC


def recallFraction(truth_h5, measured_h5, tolerance):
    """
    Return the fraction of truth localizations that have a
    measured localization within X pixels (in the XY plane).

    truth_h5 - A saH5Py.SAH5Py object with the ground truth localizations.
    measured_h5 - A saH5Py.SAH5Py object with the found localizations.
    tolerance - The search radius in pixels.
    """
    if (measured_h5.getNLocalizations() == 0):
        return [0, truth_h5.getNLocalizations()]
    
    recalled_locs = 0
    total_locs = 0
    for i in range(truth_h5.getMovieLength()):
        t_locs = truth_h5.getLocalizationsInFrame(i)
        m_locs = measured_h5.getLocalizationsInFrame(i)

        if bool(t_locs) and bool(m_locs):
            dist = iaUtilsC.peakToPeakDistAndIndex(t_locs['x'], t_locs['y'],
                                                   m_locs['x'], m_locs['y'],
                                                   max_distance = tolerance)[0]

            recalled_locs += numpy.count_nonzero((dist >= 0.0))
            total_locs += dist.size
        elif bool(t_locs):
            total_locs += t_locs['x'].size

    return [recalled_locs, total_locs]


def noiseFraction(truth_h5, measured_h5, tolerance):
    """
    Return the fraction of measured localizations that are greater than
    tolerance pixels from the nearest truth localization.

    Note: This will return 0 if there are no measured localizations.

    truth_h5 - A saH5Py.SAH5Py object with the ground truth localizations.
    measured_h5 - A saH5Py.SAH5Py object with the found localizations.
    tolerance - The search radius in pixels.
    """
    if (measured_h5.getNLocalizations() == 0):
        return [0, truth_h5.getNLocalizations()]
    
    noise_locs = 0
    total_locs = 0
    for i in range(truth_h5.getMovieLength()):
        t_locs = truth_h5.getLocalizationsInFrame(i)
        m_locs = measured_h5.getLocalizationsInFrame(i)

        if bool(t_locs) and bool(m_locs):
            dist = iaUtilsC.peakToPeakDistAndIndex(m_locs['x'], m_locs['y'],
                                                   t_locs['x'], t_locs['y'],
                                                   max_distance = tolerance)[0]

            noise_locs += numpy.count_nonzero((dist < 0.0))
            total_locs += dist.size
        elif bool(m_locs):
            noise_locs += m_locs['x'].size
            total_locs += m_locs['x'].size

    return [noise_locs, total_locs]


if (__name__ == "__main__"):

    import argparse

    import storm_analysis.sa_library.sa_h5py as saH5Py


    parser = argparse.ArgumentParser(description = 'Calculate recall fraction (in XY).')

    parser.add_argument('--truth_bin', dest='truth_bin', type=str, required=True,
                        help = "Ground truth localization file.")
    parser.add_argument('--measured_bin', dest='measured_bin', type=str, required=True,
                        help = "Measured localization file.")
    parser.add_argument('--tolerance', dest='tolerance', type=float, default = 0.2, required=False,
                        help = "Tolerance in position difference in pixels.")

    args = parser.parse_args()

    truth_h5 = saH5Py.SAH5Py(args.truth_bin)
    measured_h5 = saH5Py.SAH5Py(args.measured_bin)

    [recalled_locs, total_locs] = recallFraction(truth_h5, measured_h5, args.tolerance)
    print("Recall fraction {0:.5f}".format(float(recalled_locs)/float(total_locs)))

    [noise_locs, total_locs] = noiseFraction(truth_h5, measured_h5, args.tolerance)
    print("Noise fraction {0:.5f}".format(float(noise_locs)/float(total_locs)))

    truth_h5.close()
    measured_h5.close()

#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
