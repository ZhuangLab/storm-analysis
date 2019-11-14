#!/usr/bin/env python
"""
Perform spline analysis on a dax file given parameters.

Hazen 01/16
"""

import storm_analysis.spliner.find_peaks_decon as findPeaksDecon
import storm_analysis.spliner.find_peaks_std as findPeaksSTD

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as std_analysis


class SplineAnalysisException(Exception):
    pass


def analyze(movie_name, mlist_name, settings_name):

    # Load parameters.
    parameters = params.ParametersSpliner().initFromFile(settings_name, warnings = False)

    # Check for v1.0 parameters.
    if not (parameters.hasAttr("camera_gain") or parameters.hasAttr("camera_calibration")):
        raise SplineAnalysisException("Camera parameters are missing. Version 1.0 parameters?")
    
    # Create appropriate finding and fitting object.
    if parameters.hasAttr("decon_method"):
        finder = findPeaksDecon.initFindAndFit(parameters, settings_name)
    else:
        parameters = params.ParametersSplinerSTD().initFromFile(settings_name)
        finder = findPeaksSTD.initFindAndFit(parameters)

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
    data_writer = analysisIO.DataWriterHDF5(data_file = mlist_name,
                                            parameters = parameters,
                                            sa_type = "Spliner")
        
    std_analysis.standardAnalysis(finder,
                                  movie_reader,
                                  data_writer,
                                  parameters)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'C-Spline analysis - Babcock and Zhuang, Scientific Reports, 2017')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file. This is a binary file in Insight3 format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    analyze(args.movie, args.mlist, args.settings)

