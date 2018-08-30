#!/usr/bin/env python
"""
Estimates the static background in a STORM movie.

The estimate is performed by averaging this might
not be the best choice for movies with a high density
of real localizations.

This may be a good choice if you have a largish
fixed background and a relatively low density of
real localizations.

Hazen 8/16
"""

import numpy


class StaticBackgroundException(Exception):
    
    def __init__(self, message):
        Exception.__init__(self, message)


class StaticBGEstimator(object):
    """
    Estimates the background using a simple boxcar average.

    In the case of movies with activation frames, these frames will
    be ignored in the background estimate.

    Note: This expects to be asked for estimates in a sequential 
    fashion as would occur during normal STORM movie analysis.
    """
    def __init__(self, frame_reader = None, start_frame = 0, sample_size = 100, descriptor = "1", **kwds):
        self.cur_frame = start_frame - 1
        self.descriptor = descriptor
        self.descriptor_len = len(descriptor)
        self.frame_reader = frame_reader
        self.number_averaged = 0
        self.sample_size = sample_size

        [movie_w, movie_h, self.movie_l] = frame_reader.filmSize()

        # Figure out where to start and end the average.
        end_frame = start_frame + int(self.sample_size/2)
        start_frame = start_frame - int(self.sample_size/2)
        
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
            if not self.shouldIgnore(i):
                self.number_averaged += 1
                self.running_sum += self.frame_reader.loadAFrame(i)

    def estimateBG(self, frame_number):
        if (frame_number != (self.cur_frame + 1)):
            raise StaticBackgroundException("Received request for an estimate of a non-sequential frame " + str(self.cur_frame) + " " + str(frame_number))
        else:
            self.cur_frame = frame_number

        # Move average forward by 1 frame if possible.
        start_frame = frame_number - int(self.sample_size/2)
        end_frame = frame_number + int(self.sample_size/2)
        if (start_frame > 0) and (end_frame < self.movie_l):

            # Remove old frame.
            if not self.shouldIgnore(start_frame - 1):
                self.number_averaged -= 1
                self.running_sum -= self.frame_reader.loadAFrame(start_frame - 1)

            # Add new frame.
            if not self.shouldIgnore(end_frame):
                self.number_averaged += 1
                self.running_sum += self.frame_reader.loadAFrame(end_frame)

        # Return the current average.
        return self.running_sum/self.number_averaged

    def shouldIgnore(self, frame_number):
        desc = self.descriptor[frame_number % self.descriptor_len]
        if (desc == "0"):
            #print("Ignoring frame", frame_number)
            return True
        else:
            return False


if (__name__ == "__main__"):

    import argparse

    import storm_analysis.sa_library.datareader as datareader
    import storm_analysis.sa_library.datawriter as datawriter
    import storm_analysis.sa_library.parameters as params

    # Process command line arguments.
    parser = argparse.ArgumentParser(description = 'Running average background subtraction')

    parser.add_argument('--in_movie', dest='in_movie', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--out_movie', dest='out_movie', type=str, required=True,
                        help = "The name of the output movie (with background subtracted). This will be in .dax format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()

    # Load movies and parameters.
    input_movie = datareader.inferReader(args.in_movie)
    [w, h, l] = input_movie.filmSize()
    
    output_movie = datawriter.DaxWriter(args.out_movie)
    parameters = params.ParametersCommon().initFromFile(args.settings)
    
    n_frames = parameters.getAttr("max_frame")
    if (n_frames > l) or (n_frames == -1):
        n_frames = l

    # Default to a sample size if the settings file does not specify this.
    sample_size = 100
    if (parameters.getAttr("static_background_estimate", 0) > 0):
        sample_size = parameters.getAttr("static_background_estimate")
    else:
        print("Did not find parameter 'static_background_estimate' in parameters file, defaulting to", sample_size)
            
    sbge = StaticBGEstimator(input_movie,
                             sample_size = sample_size,
                             descriptor = parameters.getAttr("descriptor"))
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
