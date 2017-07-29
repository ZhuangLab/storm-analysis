#!/usr/bin/env python
"""
Perform multi-plane analysis on a SMLM movie given parameters.

This can also be used to perform spline based fitting to single
plane sCMOS camera data. Just comment out the 'mapping' parameter,
or set it to a file that does not exist.

Hazen 05/17
"""

import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as stdAnalysis

import storm_analysis.multi_plane.find_peaks_std as findPeaksStd


def analyze(base_name, mlist_name, settings_name):
    parameters = params.ParametersMultiplane().initFromFile(settings_name)
    finder = findPeaksStd.MPPeakFinder(parameters)
    fitter = findPeaksStd.MPPeakFitter(parameters)
    finder_fitter = findPeaksStd.MPFinderFitter(parameters, finder, fitter)
    reader = findPeaksStd.MPMovieReader(base_name = base_name,
                                        parameters = parameters)
    data_writer = findPeaksStd.MPDataWriter(data_file = mlist_name,
                                            parameters = parameters)
    stdAnalysis.standardAnalysis(finder_fitter,
                                 reader,
                                 data_writer,
                                 parameters)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Multi-plane analysis - ...')

    parser.add_argument('--basename', dest='basename', type=str, required=True,
                        help = "The base name of the movie to analyze. Movies can be in .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file. This is a binary file in Insight3 format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    analyze(args.basename, args.mlist, args.settings)
