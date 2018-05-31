#!/usr/bin/env python
"""
Reslices a camera calibration file. This is used to
create a calibration that matches the camera ROI of
the acquired data.

Hazen 05/18
"""
import numpy

import storm_analysis.sa_library.analysis_io as analysisIO


def resliceCalibration(in_cal, out_cal, xs, ys, xw, yw):
    """
    This follows the same X/Y convention as the image mask and the
    localizations, Y is the slow axis and X is the fast axis.
    """

    # Load the data & reshape.
    [offset, variance, gain, rqe] = analysisIO.loadCMOSCalibration(in_cal)

    # Slice out the ROI.
    xe = xs + xw
    ye = ys + yw
    rs_offset = offset[ys:ye,xs:xe]
    rs_variance = variance[ys:ye,xs:xe]
    rs_gain = gain[ys:ye,xs:xe]
    rs_rqe = rqe[ys:ye,xs:xe]

    # Save sliced calibration.
    numpy.save(out_cal, [rs_offset, rs_variance, rs_gain, rs_rqe, 2])


if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = 'sCMOS calibration re-slicer')

    parser.add_argument('--in_cal', dest='in_cal', type=str, required=True,
                        help = "Input calibration file.")
    parser.add_argument('--out_cal', dest='out_cal', type=str, required=True,
                        help = "Output calibration file.")
    parser.add_argument('--xs', dest='xs', type=int, required=True,
                        help = "x start (pixels).")
    parser.add_argument('--xw', dest='xw', type=int, required=True,
                        help = "x width (pixels).")
    parser.add_argument('--ys', dest='ys', type=int, required=True,
                        help = "y start (pixels).")
    parser.add_argument('--yw', dest='yw', type=int, required=True,
                        help = "y width (pixels).")

    args = parser.parse_args()    
    
    resliceCalibration(args.in_cal, args.out_cal, args.xs, args.ys, args.xw, args.yw)


#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
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
