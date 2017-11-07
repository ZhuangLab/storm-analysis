#!/usr/bin/env python
"""
Functions for calculating recall, etc.

Hazen 09/17
"""
import numpy

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC


def recallFraction(truth_i3, measured_i3, tolerance):
    """
    Return the fraction of truth localizations that have a
    measured localization within X pixels (in the XY plane).

    truth_i3 - A readinsight3.I3Reader object with the ground truth localizations.
    measured_i3 - A readinsight3.I3Reader object with the found localizations.
    tolerance - The search radius in pixels.
    """
    if (measured_i3.getNumberMolecules() == 0):
        return [0, truth_i3.getNumberMolecules()]
    
    recalled_locs = 0
    total_locs = 0
    for i in range(truth_i3.getNumberFrames()):
        t_locs = truth_i3.getMoleculesInFrame(i+1)
        m_locs = measured_i3.getMoleculesInFrame(i+1, good_only = False)
        
        dist = iaUtilsC.peakToPeakDistAndIndex(t_locs['xc'], t_locs['yc'],
                                               m_locs['xc'], m_locs['yc'])[0]

        recalled_locs += numpy.count_nonzero((dist < tolerance))
        total_locs += dist.size

    return [recalled_locs, total_locs]


def noiseFraction(truth_i3, measured_i3, tolerance):
    """
    Return the fraction of measured localizations that are greater than
    tolerance pixels from the nearest truth localization.

    Note: This will return 0 if there are no measured localizations.

    truth_i3 - A readinsight3.I3Reader object with the ground truth localizations.
    measured_i3 - A readinsight3.I3Reader object with the found localizations.
    tolerance - The search radius in pixels.
    """
    if (measured_i3.getNumberMolecules() == 0):
        return [0, truth_i3.getNumberMolecules()]
    
    noise_locs = 0
    total_locs = 0
    for i in range(truth_i3.getNumberFrames()):
        t_locs = truth_i3.getMoleculesInFrame(i+1)
        m_locs = measured_i3.getMoleculesInFrame(i+1, good_only = False)

        dist = iaUtilsC.peakToPeakDistAndIndex(t_locs['xc'], t_locs['yc'],
                                               m_locs['xc'], m_locs['yc'])[0]

        noise_locs += numpy.count_nonzero((dist > tolerance))
        total_locs += dist.size

    return [noise_locs, total_locs]


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Calculate recall fraction (in XY).')

    parser.add_argument('--truth_bin', dest='truth_bin', type=str, required=True,
                        help = "Ground truth localization file.")
    parser.add_argument('--measured_bin', dest='measured_bin', type=str, required=True,
                        help = "Measured localization file.")
    parser.add_argument('--tolerance', dest='tolerance', type=float, default = 0.2, required=False,
                        help = "Tolerance in position difference in pixels.")

    args = parser.parse_args()

    truth_i3 = readinsight3.I3Reader(args.truth_bin)
    measured_i3 = readinsight3.I3Reader(args.measured_bin)

    [recalled_locs, total_locs] = recallFraction(truth_i3, measured_i3, args.tolerance)
    print("Recall fraction {0:.5f}".format(float(recalled_locs)/float(total_locs)))

    [noise_locs, total_locs] = noiseFraction(truth_i3, measured_i3, args.tolerance)
    print("Noise fraction {0:.5f}".format(float(noise_locs)/float(total_locs)))

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
