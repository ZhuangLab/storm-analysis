#!/usr/bin/env python
"""
Perform sCMOS analysis on a SMLM movie given parameters.

Note: This now differs somewhat from the original sCMOS algorithm by Huang 
      in that it no longer does the peak finding by convolution with a zero-
      mean PSF. Instead it follows the approach described in this paper:

      "SNSMIL, a real-time single molecule identification and localization 
      algorithm for super-resolution fluorescence microscopy", Tang et al., 
      Scientific Reports, 2015.

Hazen 09/17
"""

import storm_analysis.daostorm_3d.find_peaks as find_peaks

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as std_analysis


def analyze(movie_name, mlist_name, settings_name):

    # Load parameters.
    parameters = params.ParametersSCMOS().initFromFile(settings_name)

    # Check for possibly v1.0 parameters.
    if not parameters.hasAttr("background_sigma"):
        raise Exception("Parameter 'background_sigma' is missing. Version 1.0 parameters?")
    
    # Create finding and fitting object.
    finder = find_peaks.initFindAndFit(parameters)

    # Create object for reading sCMOS camera frames.
    frame_reader = analysisIO.FrameReaderSCMOS(movie_file = movie_name,
                                               parameters = parameters)

    # Create movie reader (uses frame reader).
    movie_reader = analysisIO.MovieReader(frame_reader = frame_reader,
                                          parameters = parameters)
    
    # Create localization file writer.
    data_writer = analysisIO.DataWriterHDF5(data_file = mlist_name,
                                            parameters = parameters,
                                            sa_type = 'sCMOS')

    # Run the analysis.
    std_analysis.standardAnalysis(finder,
                                  movie_reader,
                                  data_writer,
                                  parameters)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'sCMOS analysis following Huang, Nature Methods, 2013')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file. This is a binary file in Insight3 format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    analyze(args.movie, args.mlist, args.settings)


#
# The MIT License
#
# Copyright (c) 2014 Zhuang Lab, Harvard University
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
