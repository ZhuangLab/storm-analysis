#!/usr/bin/env python
"""
Perform multi-plane 3D-DAOSTORM analysis on a SMLM movie given parameters.

Hazen 01/18
"""

import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as stdAnalysis

import storm_analysis.multi_plane.analysis_io as analysisIO
import storm_analysis.multi_plane.find_peaks_dao_std as findPeaksStd


def analyze(base_name, mlist_name, settings_name):

    # Load parameters.
    parameters = params.ParametersMultiplaneDao().initFromFile(settings_name)

    # Create finding and fitting object.
    finder = findPeaksStd.initFindAndFit(parameters)

    # Create multiplane (sCMOS) reader.
    reader = analysisIO.MPMovieReader(base_name = base_name,
                                      parameters = parameters)

    # Create multiplane localization file(s) writer.
    data_writer = analysisIO.MPDataWriter(data_file = mlist_name,
                                          parameters = parameters,
                                          sa_type = "Multiplane-Dao")

    stdAnalysis.standardAnalysis(finder,
                                 reader,
                                 data_writer,
                                 parameters)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Multi-plane 3D-DAOSTORM analysis - ...')

    parser.add_argument('--basename', dest='basename', type=str, required=True,
                        help = "The base name of the movie to analyze. Movies can be in .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    analyze(args.basename, args.mlist, args.settings)
