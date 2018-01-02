#!/usr/bin/env python
"""
Performs "standard" analysis on a dax file given parameters.

Hazen 10/13
"""

import numpy
import os

import storm_analysis.sa_utilities.apply_drift_correction_c as applyDriftCorrectionC
import storm_analysis.sa_utilities.avemlist_c as avemlistC
import storm_analysis.sa_utilities.fitz_c as fitzC
import storm_analysis.sa_utilities.tracker as tracker
import storm_analysis.sa_utilities.xyz_drift_correction as xyzDriftCorrection

                    
def averaging(mol_list_filename, ave_list_filename):
    """
    Averages all the molecules in a track into a single molecule.
    """
    avemlistC.avemlist(mol_list_filename, ave_list_filename)

def driftCorrection(list_files, parameters):
    """
    Performs drift correction.
    """
    drift_name = list_files[0][:-9] + "drift.txt"

    # Check if we have been asked not to do z drift correction.
    # The default is to do the correction.
    z_correct = True
    if (parameters.getAttr("z_correction", 0) != 0):
        z_correct = False

    #
    # Get z range from the paraemeters file. Note these are
    # in microns and we are using nanometers.
    #
    [min_z, max_z] = parameters.getZRange()
            
    xyzDriftCorrection.xyzDriftCorrection(list_files[0],
                                          drift_name,
                                          parameters.getAttr("frame_step"),
                                          parameters.getAttr("d_scale"),
                                          1000.0 * min_z,
                                          1000.0 * max_z,
                                          z_correct)

    if (os.path.exists(drift_name)):
        for list_file in list_files:
            applyDriftCorrectionC.applyDriftCorrection(list_file, drift_name)

def peakFinding(find_peaks, movie_reader, data_writer, parameters):
    """
    Does the peak finding.
    """
    curf = data_writer.getStartFrame()
    movie_reader.setup(curf)

    #
    # Analyze the movie.
    #
    # Catch keyboard interrupts & "gracefully" exit.
    #
    try:
        while(movie_reader.nextFrame()):

            # Find the localizations.
            peaks = find_peaks.analyzeImage(movie_reader)

            # Save results
            data_writer.addPeaks(peaks, movie_reader)

            print("Frame:",
                  movie_reader.getCurrentFrameNumber(),
                  data_writer.getNumberAdded(),
                  data_writer.getTotalPeaks())

        print("")
        data_writer.close()
        find_peaks.cleanUp()
        return True

    except KeyboardInterrupt:
        print("Analysis stopped.")
        data_writer.close()
        find_peaks.cleanUp()
        return False

def standardAnalysis(find_peaks, movie_reader, data_writer, parameters):
    """
    Perform standard analysis.

    movie_reader - sa_utilities.analysis_io.MovieReader object.
    data_writer - sa_utilities.analysis_io.DataWriter object.
    """
    # Peak finding
    print("Peak finding")
    if peakFinding(find_peaks, movie_reader, data_writer, parameters):

        # Z fitting, '3d' model, localizations.
        if parameters.getAttr("do_zfit", False):
            if (parameters.getAttr("model", "") == "3d"):
                print("Localization z fitting")
                zFitting(data_writer.getFilename(), parameters, False)
                
        # Drift correction.
        print("")
        mlist_file = data_writer.getFilename()

        # Tracking and averaging.
        print("Tracking")
        tracker.tracker(data_writer.getFilename(),
                        descriptor = parameters.getAttr("descriptor"),
                        max_gap = parameters.getAttr("max_gap", 0),
                        radius = parameters.getAttr("radius"))

        # Z fitting, '3d' model, tracks.
        if parameters.getAttr("do_zfit", False):
            if (parameters.getAttr("model", "") == "3d"):
                print("Tracks z fitting")
                zFitting(data_writer.getFilename(), parameters, True)

        # Mark out of z range tracks as category 9.

    print("Analysis complete")

def zFitting(h5_name, parameters, fit_tracks):
    """
    Does z fitting for the '3d' model.
    """
    [wx_params, wy_params] = parameters.getWidthParams()
    [min_z, max_z] = parameters.getZRange()
    z_step = parameters.getAttr("z_step", 1.0e-3)

    if fit_tracks:
        fitzC.fitzTracks(h5_name,
                         parameters.getAttr("cutoff"),
                         wx_params,
                         wy_params,
                         min_z,
                         max_z,
                         z_step)
    else:
        fitzC.fitzRaw(h5_name,
                      parameters.getAttr("cutoff"),
                      wx_params,
                      wy_params,
                      min_z,
                      max_z,
                      z_step)


#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
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
