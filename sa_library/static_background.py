#!/usr/bin/env python
#
# Estimates the static background in a STORM movie.
#
# The estimate is performed by averaging this might
# not be the best choice for movies with a high density
# of real localizations.
#
# This may be a good choice if you have a largish
# fixed background and a relatively low density of
# real localizations.
#
# Hazen 8/16
#

import numpy


class StaticBackgroundException(Exception):
    
    def __init__(self, message):
        Exception.__init__(self, message)


class StaticBGEstimator(object):
    """
    Estimates the background using a simple boxcar average. 

    Note: This expects to be asked for estimates in a sequential 
    fashion as would occur during normal STORM movie analysis.
    """
    def __init__(self, movie_data, start_frame = 0, sample_size = 100):
        self.cur_frame = start_frame - 1
        self.movie_data = movie_data
        self.sample_size = sample_size

        [movie_w, movie_h, self.movie_l] = movie_data.filmSize()

        # Figure out where to start and end the average.
        end_frame = start_frame + self.sample_size/2
        start_frame = start_frame - self.sample_size/2
        
        if (start_frame < 0):
            start_frame = 0
            end_frame = start_frame + self.sample_size
            if (end_frame > self.movie_l):
                end_frame = self.movie_l
                self.sample_size = self.movie_l

        if (end_frame > self.movie_l):
            end_frame = self.movie_l
            start_frame = end_frame - self.sample_size
            if (start_frame < 0):
                start_frame = 0
                self.sample_size = self.movie_l

        self.running_sum = numpy.zeros((movie_h, movie_w))
        for i in range(start_frame, end_frame):
            self.running_sum += self.movie_data.loadAFrame(i)

    def estimateBG(self, frame_number):
        if (frame_number != (self.cur_frame + 1)):
            raise StaticBackgroundException("Received request for an estimate of a non-sequential frame " + str(self.cur_frame) + " " + str(frame_number))
        else:
            self.cur_frame = frame_number

        # Move average forward by 1 frame if possible.
        start_frame = frame_number - self.sample_size/2
        end_frame = frame_number + self.sample_size/2
        if (start_frame > 0) and (end_frame < self.movie_l):

            # Remove old frame.
            self.running_sum -= self.movie_data.loadAFrame(start_frame - 1)

            # Add new frame.
            self.running_sum += self.movie_data.loadAFrame(end_frame)

        # Return the current average.
        return self.running_sum/self.sample_size


if (__name__ == "__main__"):

    import sys

    import sa_library.datareader as datareader
    import sa_library.daxwriter as daxwriter

    if (len(sys.argv) != 4):
        print "usage: <input movie> <output movie> <number of frames>"
        exit()

    input_movie = datareader.inferReader(sys.argv[1])
    [w, h, l] = input_movie.filmSize()
    
    output_movie = daxwriter.DaxWriter(sys.argv[2], w, h)

    n_frames = int(sys.argv[3])
    if (n_frames > l):
        n_frames = l

    sbge = StaticBGEstimator(input_movie)
    for i in range(n_frames):
        diff = input_movie.loadAFrame(i) - sbge.estimateBG(i) + 100
        output_movie.addFrame(diff)

    output_movie.close()
    
    
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
