#!/usr/bin/env python
"""
Return the fraction of truth localizations that have a
measured localization within X pixels (in the XY plane).

Hazen 04/16
"""
import numpy

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.ia_utilities_c as utilC


def recallFraction(truth_name, measured_name, tolerance):
    truth_i3 = readinsight3.I3Reader(truth_name)
    measured_i3 = readinsight3.I3Reader(measured_name)

    recalled_locs = 0
    total_locs = 0
    for i in range(truth_i3.getNumberFrames()):
        t_locs = truth_i3.getMoleculesInFrame(i+1)
        m_locs = measured_i3.getMoleculesInFrame(i+1, good_only = False)
        
        dist = utilC.peakToPeakDist(t_locs['xc'], t_locs['yc'], m_locs['xc'], m_locs['yc'])

        recalled_locs += numpy.count_nonzero((dist < tolerance))
        total_locs += dist.size

    return [recalled_locs, total_locs]


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

    [recalled_locs, total_locs] = recallFraction(args.truth_bin, args.measured_bin, args.tolerance)

    print("Recall fraction {0:.5f}".format(float(recalled_locs)/float(total_locs)))


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
