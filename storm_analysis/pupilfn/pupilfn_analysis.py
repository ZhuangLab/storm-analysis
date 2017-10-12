#!/usr/bin/env python
"""
Perform pupil function analysis on a SMLM movie given parameters.

Hazen 10/17
"""

import storm_analysis.pupilfn.find_peaks_std as find_peaks_std

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as std_analysis


def analyze(movie_name, mlist_name, settings_name):

    # Load parameters.
    parameters = params.ParametersPupilFn().initFromFile(settings_name, warnings = False)

    # Create finding and fitting object.
    finder = find_peaks_std.initFindAndFit(parameters)

    # Create appropriate reader.
    if parameters.hasAttr("camera_offset"):
        frame_reader = analysisIO.FrameReaderStd(movie_file = movie_name,
                                                 parameters = parameters)
    else:
        frame_reader = analysisIO.FrameReaderSCMOS(movie_file = movie_name,
                                                   parameters = parameters)

    # Create movie reader (uses frame reader).
    movie_reader = analysisIO.MovieReader(frame_reader = frame_reader,
                                          parameters = parameters)
    
    # Create localization file writer.
    data_writer = analysisIO.DataWriter(data_file = mlist_name,
                                        parameters = parameters)
        
    std_analysis.standardAnalysis(finder,
                                  movie_reader,
                                  data_writer,
                                  parameters)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Pupil function analysis - ...')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file. This is a binary file in Insight3 format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    analyze(args.movie, args.mlist, args.settings)

