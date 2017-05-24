#!/usr/bin/env python
"""
Perform spline analysis on a dax file given parameters.

Hazen 01/16
"""

import storm_analysis.spliner.find_peaks_fista as find_peaks_fista
import storm_analysis.spliner.find_peaks_std as find_peaks_std
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as std_analysis


def analyze(movie_name, mlist_name, settings_name):
    parameters = params.ParametersSpliner().initFromFile(settings_name, warnings = False)
    if (parameters.getAttr("use_fista", 0) != 0):
        parameters = params.ParametersSplinerFISTA().initFromFile(settings_name)
        finder = find_peaks_fista.SplinerFISTAFinderFitter(parameters)
    else:
        parameters = params.ParametersSplinerSTD().initFromFile(settings_name)
        finder = find_peaks_std.SplinerFinderFitter(parameters)        
    std_analysis.standardAnalysis(finder,
                                  movie_name,
                                  mlist_name,
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

