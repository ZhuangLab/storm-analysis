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
import storm_analysis.sa_utilities.tracker_c as trackerC
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
    # peak finding
    print("Peak finding")
    peakFinding(find_peaks, movie_reader, data_writer, parameters)
            
    if False:
        print("")
        mlist_file = data_writer.getFilename()
        
        # tracking
        print("Tracking")
        tracking(mlist_file, parameters)

        # averaging
        alist_file = None
        if (parameters.getAttr("radius") > 0.0):
            alist_file = mlist_file[:-9] + "alist.bin"
            averaging(mlist_file, alist_file)
            print("")

        # z fitting
        if (parameters.getAttr("do_zfit", 0) != 0):
            if (parameters.getAttr("model", "") == "3d"):
                print("Fitting Z")
                if alist_file:
                    zFitting(alist_file, parameters)
                zFitting(mlist_file, parameters)
                print("")
            else:
                print("Warning! Ignoring 'do_zfit' because fitting model is not '3d'!")
                print("")

        # drift correction
        if (parameters.getAttr("drift_correction", 0) != 0):
            print("Drift Correction")
            if alist_file:
                driftCorrection([mlist_file, alist_file], parameters)
            else:
                driftCorrection([mlist_file], parameters)
            print("")
    print("Analysis complete")

def tracking(mol_list_filename, parameters):
    """
    Does the frame-to-frame tracking.
    """
    [min_z, max_z] = parameters.getZRange()
    trackerC.tracker(mol_list_filename,
                     parameters.getAttr("descriptor"),
                     parameters.getAttr("radius"),
                     1000.0*min_z, 1000.0*max_z, 1)

def zFitting(mol_list_filename, parameters):
    """
    Does z fitting.
    """
    [wx_params, wy_params] = parameters.getWidthParams()
    [min_z, max_z] = parameters.getZRange()
    fitzC.fitz(mol_list_filename,
               parameters.getAttr("cutoff"),
               wx_params,
               wy_params,
               min_z * 1000.0,
               max_z * 1000.0,
               parameters.getAttr("z_step", 1.0))


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
