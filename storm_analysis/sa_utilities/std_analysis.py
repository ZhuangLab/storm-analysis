#!/usr/bin/python
#
# Performs "standard" analysis on a dax file given parameters.
#
# Hazen 10/13
#

import numpy
import os
import subprocess

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.static_background as static_background
import storm_analysis.sa_library.writeinsight3 as writeinsight3

import storm_analysis.sa_utilities.apply_drift_correction_c as applyDriftCorrectionC
import storm_analysis.sa_utilities.avemlist_c as avemlistC
import storm_analysis.sa_utilities.fitz_c as fitzC
import storm_analysis.sa_utilities.tracker_c as trackerC
import storm_analysis.sa_utilities.xyz_drift_correction as xyzDriftCorrection


src_dir = os.path.dirname(__file__)
if not (src_dir == ""):
    src_dir += "/"

# Averages all the molecules in a track into a single molecule.
def averaging(mol_list_filename, ave_list_filename):
    avemlistC.avemlist(mol_list_filename, ave_list_filename)

# Performs drift correction.
def driftCorrection(list_files, parameters):
    drift_name = list_files[0][:-9] + "drift.txt"

    # Check if we have been asked not to do z drift correction.
    # The default is to do the correction.
    z_correct = True
    if hasattr(parameters, "z_correction"):
        if not parameters.z_correction:
            z_correct = False

    xyzDriftCorrection.xyzDriftCorrection(list_files[0],
                                          drift_name,
                                          parameters.frame_step,
                                          parameters.d_scale,
                                          correct_z = z_correct)

    if (os.path.exists(drift_name)):
        for list_file in list_files:
            applyDriftCorrectionC.applyDriftCorrection(list_file, drift_name)

# Does the peak finding.
def peakFinding(find_peaks, movie_file, mlist_file, parameters):

    # open files for input & output
    movie_data = datareader.inferReader(movie_file)
    [movie_x,movie_y,movie_l] = movie_data.filmSize()

    # if the i3 file already exists, read it in,
    # write it out & start the analysis from the
    # end.
    total_peaks = 0
    if(os.path.exists(mlist_file)):
        print("Found", mlist_file)
        i3data_in = readinsight3.loadI3File(mlist_file)
        try:
            curf = int(numpy.max(i3data_in['fr']))
        except ValueError:
            curf = 0
        print(" Starting analysis at frame:", curf)
        i3data = writeinsight3.I3Writer(mlist_file)
        if (curf > 0):
            i3data.addMolecules(i3data_in)
            total_peaks = i3data_in['x'].size
    else:
        curf = 0
        i3data = writeinsight3.I3Writer(mlist_file)

    # process parameters
    if hasattr(parameters, "start_frame"):
        if (parameters.start_frame>=curf) and (parameters.start_frame<movie_l):
            curf = parameters.start_frame

    if hasattr(parameters, "max_frame"):
        if (parameters.max_frame>0) and (parameters.max_frame<movie_l):
            movie_l = parameters.max_frame

    static_bg_estimator = None
    if hasattr(parameters, "static_background_estimate"):
        if (parameters.static_background_estimate > 0):
            print("Using static background estimator.")
            static_bg_estimator = static_background.StaticBGEstimator(movie_data,
                                                                      start_frame = curf,
                                                                      sample_size = parameters.static_background_estimate)

    # analyze the movie
    # catch keyboard interrupts & "gracefully" exit.
    try:
        while(curf<movie_l):
            #for j in range(l):

            # Set up the analysis.
            image = movie_data.loadAFrame(curf) - parameters.baseline
            mask = (image < 1.0)
            if (numpy.sum(mask) > 0):
                print(" Removing negative values in frame", curf)
                image[mask] = 1.0

            # Find and fit the peaks.
            if static_bg_estimator is not None:
                bg_estimate = static_bg_estimator.estimateBG(curf) - parameters.baseline
                [peaks, residual] = find_peaks.analyzeImage(image,
                                                            bg_estimate = bg_estimate)
            else:
                [peaks, residual] = find_peaks.analyzeImage(image)

            # Save the peaks.
            if (type(peaks) == type(numpy.array([]))):
                # remove unconverged peaks
                peaks = find_peaks.getConvergedPeaks(peaks)

                # save results
                if(parameters.orientation == "inverted"):
                    i3data.addMultiFitMolecules(peaks, movie_x, movie_y, curf+1, parameters.pixel_size, inverted = True)
                else:
                    i3data.addMultiFitMolecules(peaks, movie_x, movie_y, curf+1, parameters.pixel_size, inverted = False)

                total_peaks += peaks.shape[0]
                print("Frame:", curf, peaks.shape[0], total_peaks)
            else:
                print("Frame:", curf, 0, total_peaks)
            curf += 1

        print("")
        i3data.close()
        find_peaks.cleanUp()
        return 0

    except KeyboardInterrupt:
        print("Analysis stopped.")
        i3data.close()
        find_peaks.cleanUp()
        return 1

# Perform standard analysis.
def standardAnalysis(find_peaks, data_file, mlist_file, parameters):

    # peak finding
    print("Peak finding")
    if(not peakFinding(find_peaks, data_file, mlist_file, parameters)):
        print("")
        
        # tracking
        print("Tracking")
        tracking(mlist_file, parameters)

        # averaging
        alist_file = None
        if(parameters.radius > 0.0):
            alist_file = mlist_file[:-9] + "alist.bin"
            averaging(mlist_file, alist_file)
            print("")

        # z fitting
        if hasattr(parameters, "do_zfit") and parameters.do_zfit:
            print("Fitting Z")
            if alist_file:
                zFitting(alist_file, parameters)
            zFitting(mlist_file, parameters)
            print("")

        # drift correction
        if hasattr(parameters, "drift_correction"):
            if parameters.drift_correction:
                print("Drift Correction")
                if alist_file:
                    driftCorrection([mlist_file, alist_file], parameters)
                else:
                    driftCorrection([mlist_file], parameters)
                print("")
    print("Analysis complete")

# Does the frame-to-frame tracking.
def tracking(mol_list_filename, parameters):
    [min_z, max_z] = params.getZRange(parameters)
    trackerC.tracker(mol_list_filename, parameters.descriptor, parameters.radius, 1000.0*min_z, 1000.0*max_z, 1)

# Does z fitting.
def zFitting(mol_list_filename, parameters):
    if(parameters.orientation == "inverted"):
        wx_params = params.getWidthParams(parameters, "y")
        wy_params = params.getWidthParams(parameters, "x")
    else:
        wx_params = params.getWidthParams(parameters, "x")
        wy_params = params.getWidthParams(parameters, "y")
    fitzC.fitz(mol_list_filename, parameters.cutoff, wx_params, wy_params)

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
