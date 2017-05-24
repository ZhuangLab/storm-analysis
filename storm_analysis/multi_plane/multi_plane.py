#!/usr/bin/env python
"""
Perform multi-plane analysis on a SMLM movie given parameters.

Hazen 05/17
"""

import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as std_analysis

import storm_analysis.multi_plane.find_peaks_std as findPeaksStd


def analyze(base_name, mlist_name, settings_name):
    parameters = params.ParametersMultiplane().initFromFile(settings_name)
    
    # And unfortunately this is different enough that  we can't just do
    # the standard analysis..


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Multi-plane analysis - ...')

    parser.add_argument('--base_name', dest='basename', type=str, required=True,
                        help = "The base name of the movie to analyze. Movie can be .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file. This is a binary file in Insight3 format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    analyze(args.movie, args.mlist, args.settings)
