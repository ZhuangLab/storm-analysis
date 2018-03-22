#!/usr/bin/env python
"""
Perform pupil function analysis on a SMLM movie given parameters.

This was written to verify our claim that a cubic spline (CS) representation of
the PSF is faster than a pupil function (PF) representation. This appears to be 
true, but not for the reasons we originally thought (that PF requires lots of
math, sines, cosines, etc..). Basically the problem with PF fitting is that the 
AOI has to be substantially larger than that of CS to achieve a similar accuracy, 
especially in Z (for astigmatism imaging). In the tested size ranges, PF was 
actually faster at any given size, but PF needs a AOI that is about 50% larger 
(linear dimension) to get the same performance. So CS at 20x20 pixels performs 
as well as PF at 30x30 pixels, and since 20x20 pixels is substantially smaller 
than 30x30 pixels CS is approximately 2x faster.

See diagnostics/spliner and diagnostics/pupilfn for more detailed results.

Currently this does not support and OTF scaling factor as described in the
Hanser paper.

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
    data_writer = analysisIO.DataWriterHDF5(data_file = mlist_name,
                                            parameters = parameters,
                                            sa_type = "Pupil-Function")

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

