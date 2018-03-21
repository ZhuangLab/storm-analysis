#!/usr/bin/env python
"""
Converts a dax file into "calibration" format, i.e. a format
that can be read by camera_calibration.py for the purpose
of camera calibration.

Note that this is more for testing. Using "calibrate.py" in
the STORM-Control project is more efficient as you don't
have to save a massive movie first.

Hazen 10/13
"""

import numpy
import sys

import storm_analysis.sa_library.datareader as datareader

def movieToCalibration(movie_name):
    """
    Calculate calibration data from a movie. This includes
    the mean intensity per frame to reduce issues with
    average brightness of the light source drifting during
    the movie.

    movie_name - The name of the movie.
    """

    # Open the movie.
    in_file = datareader.inferReader(movie_name)
    [w, h, l] = in_file.filmSize()

    # Calculate frame mean, x & xx.
    frame_mean = numpy.zeros(l)
    N = numpy.zeros((h,w), dtype = numpy.int64)
    NN = numpy.zeros((h,w), dtype = numpy.int64)

    for i in range(l):
        aframe = in_file.loadAFrame(i)
        frame_mean[i] = numpy.mean(aframe)
        
        aframe = aframe.astype(numpy.int64)
        N += aframe
        NN += aframe * aframe

    return [frame_mean, N, NN]


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Create a sCMOS calibration input file from a movie.')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie.")
    parser.add_argument('--cal', dest='cal', type=str, required=True,
                        help = "The name of the calibration input file.")

    args = parser.parse_args()

    [frame_mean, N, NN] = movieToCalibration(args.movie)

    # Save the results.
    numpy.save(args.cal, [frame_mean, N, NN])

    mean = N/float(frame_mean.size)
    print("Mean:", numpy.mean(mean))
    print("Variance:", numpy.mean(NN/float(frame_mean.size) - mean*mean))


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
