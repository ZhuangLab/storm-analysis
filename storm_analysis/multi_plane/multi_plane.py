#!/usr/bin/env python
"""
Perform multi-plane analysis on a SMLM movie given parameters.

Hazen 05/17
"""

import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as stdAnalysis

import storm_analysis.multi_plane.find_peaks_std as findPeaksStd


def analyze(base_name, mlist_name, settings_name):
    parameters = params.ParametersMultiplane().initFromFile(settings_name)
    finder = findPeaksStd.MPFinderFitter(parameters)
    reader = findPeaksStd.MPMovieReader(base_name = base_name,
                                        parameters = parameters)
    stdAnalysis.standardAnalysis(finder,
                                 reader,
                                 mlist_name,
                                 parameters)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Multi-plane analysis - ...')

    parser.add_argument('--base_name', dest='basename', type=str, required=True,
                        help = "The base name of the movie to analyze. Movies can be in .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file. This is a binary file in Insight3 format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    analyze(args.basename, args.mlist, args.settings)
