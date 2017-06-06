#!/usr/bin/env python
"""
Performs "standard" analysis on a dax file given parameters.

Hazen 10/13
"""

import numpy
import os

from xml.etree import ElementTree


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


class MovieReader(object):
    """
    Encapsulate getting the frames of a movie from a datareader.Reader, handling
    static_background, which frame to start on, etc...

    The primary reason for using a class like this is to add the flexibility
    necessary in order to make the peakFinding() function also work with
    multi-channel data / analysis.
    """
    def __init__(self, movie_file = None, parameters = None, **kwds):
        super().__init__(**kwds)

        self.background = None
        self.baseline = parameters.getAttr("baseline")
        self.bg_estimator = None
        self.cur_frame = 0
        self.frame = None
        self.max_frame = None
        self.movie_data = datareader.inferReader(movie_file)
        [self.movie_x, self.movie_y, self.movie_l] = self.movie_data.filmSize()
        self.parameters = parameters

    def getBackground(self):
        return self.background

    def getCurrentFrameNumber(self):
        return self.cur_frame
    
    def getFrame(self):
        return self.frame

    def getMovieL(self):
        return self.movie_l
    
    def getMovieX(self):
        return self.movie_x

    def getMovieY(self):
        return self.movie_y

    def hashID(self):
        return self.movie_data.hashID()

    def nextFrame(self):
        if (self.cur_frame < self.max_frame):

            # Update background estimate.
            if self.bg_estimator is not None:
                self.background = self.bg_estimator.estimateBG(self.cur_frame) - self.baseline

            # Load frame & remove all values less than 1.0 as we are doing MLE fitting.
            self.frame = self.movie_data.loadAFrame(self.cur_frame) - self.baseline
            mask = (self.frame < 1.0)
            if (numpy.sum(mask) > 0):
                print(" Removing negative value in frame", self.cur_frame)
                self.frame[mask] = 1.0

            #
            # Increment here because the .bin files are 1 indexed, but the movies
            # are 0 indexed. We'll get the right value when we call getCurrentFrameNumber().
            #
            self.cur_frame += 1
            return True
        else:
            return False
        
    def setup(self, start_frame):

        # Figure out where to start.
        self.cur_frame = start_frame
        if self.parameters.hasAttr("start_frame"):
            if (self.parameters.getAttr("start_frame") >= self.cur_frame):
                if (self.parameters.getAttr("start_frame") < self.movie_l):
                    self.cur_frame = self.parameters.getAttr("start_frame")
        
        # Figure out where to stop.
        self.max_frame = self.movie_l
        if self.parameters.hasAttr("max_frame"):
            if (self.parameters.getAttr("max_frame") > 0):
                if (self.parameters.getAttr("max_frame") < self.movie_l):
                    self.max_frame = self.parameters.getAttr("max_frame")

        # Configure background estimator, if any.
        if (self.parameters.getAttr("static_background_estimate", 0) > 0):
            print("Using static background estimator.")
            s_size = self.parameters.getAttr("static_background_estimate")
            self.bg_estimator = static_background.StaticBGEstimator(self.movie_data,
                                                                    start_frame = self.cur_frame,
                                                                    sample_size = s_size)

                    
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

    xyzDriftCorrection.xyzDriftCorrection(list_files[0],
                                          drift_name,
                                          parameters.getAttr("frame_step"),
                                          parameters.getAttr("d_scale"),
                                          correct_z = z_correct)

    if (os.path.exists(drift_name)):
        for list_file in list_files:
            applyDriftCorrectionC.applyDriftCorrection(list_file, drift_name)

def peakFinding(find_peaks, movie_reader, mlist_file, parameters):
    """
    Does the peak finding.
    """

    #
    # If the i3 file already exists, read it in, write it
    # out & start the analysis from the end.
    #
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

    movie_reader.setup(curf)

    #
    # Analyze the movie.
    #
    # Catch keyboard interrupts & "gracefully" exit.
    #
    try:
        while(movie_reader.nextFrame()):

            # Find the localizations.
            [peaks, residual] = find_peaks.analyzeImage(movie_reader)

            # Save the localizations.
            if isinstance(peaks, numpy.ndarray):
                
                # Remove unconverged localizations.
                peaks = find_peaks.getConvergedPeaks(peaks)

                # save results
                i3data.addMultiFitMolecules(peaks,
                                            movie_reader.getMovieX(),
                                            movie_reader.getMovieY(),
                                            movie_reader.getCurrentFrameNumber(),
                                            parameters.getAttr("pixel_size"),
                                            inverted = (parameters.getAttr("orientation", "normal") == "inverted"))

                total_peaks += peaks.shape[0]
                print("Frame:", movie_reader.getCurrentFrameNumber(), peaks.shape[0], total_peaks)
            else:
                print("Frame:", movie_reader.getCurrentFrameNumber(), 0, total_peaks)

        print("")
        if parameters.getAttr("append_metadata", True):

            etree = ElementTree.Element("xml")

            # Add analysis parameters.
            etree.append(parameters.toXMLElementTree())

            # Add movie properties.
            movie_props = ElementTree.SubElement(etree, "movie")
            field = ElementTree.SubElement(movie_props, "hash_value")
            field.text = movie_reader.hashID()
            for elt in [["movie_x", movie_reader.getMovieX()],
                        ["movie_y", movie_reader.getMovieY()],
                        ["movie_l", movie_reader.getMovieL()]]:
                field = ElementTree.SubElement(movie_props, elt[0])
                field.text = str(elt[1])

            i3data.closeWithMetadata(ElementTree.tostring(etree, 'ISO-8859-1'))
        else:
            i3data.close()
        find_peaks.cleanUp()
        return 0

    except KeyboardInterrupt:
        print("Analysis stopped.")
        i3data.close()
        find_peaks.cleanUp()
        return 1

def standardAnalysis(find_peaks, movie_data, mlist_file, parameters):
    """
    Perform standard analysis.

    Note: movie_data can either be a string specifying the name of a movie
          or a MovieReader object. If it is a string then a MovieReader
          object will be created.
    """
    if not isinstance(movie_data, MovieReader):
        movie_data = MovieReader(movie_file = movie_data,
                                 parameters = parameters)

    # peak finding
    print("Peak finding")
    if(not peakFinding(find_peaks, movie_data, mlist_file, parameters)):
        print("")
        
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
            print("Fitting Z")
            if alist_file:
                zFitting(alist_file, parameters)
            zFitting(mlist_file, parameters)
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
