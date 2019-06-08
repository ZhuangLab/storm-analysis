#!/usr/bin/env python
"""
Analysis IO specialized for multiplane fitting.

Hazen 09/17
"""
import numpy
import os
import sys

from xml.etree import ElementTree

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.multi_plane.mp_utilities as mpUtil


class MPDataWriter(analysisIO.DataWriter):
    """
    Data writer specialized for multi-plane data.
    """
    def __init__(self, parameters = None, sa_type = None, **kwds):
        super(MPDataWriter, self).__init__(**kwds)

        self.movie_info_set = False
        self.offsets = []

        # Figure out how many planes there are.
        self.n_planes = len(mpUtil.getExtAttrs(parameters))

        # Save frame offsets for each plane.
        for offset in mpUtil.getOffsetAttrs(parameters):
            self.offsets.append(parameters.getAttr(offset))        

        # Figure out where to start if the analysis file already exists.
        if os.path.exists(self.filename):
            print("Existing analysis file found. Restarting from last analyzed frame.")
            self.h5 = saH5Py.SAH5Py(filename = self.filename)

            self.movie_info_set = True
            
            # Find the last frame that we analyzed.
            i = self.h5.getMovieLength()
            while (i > 0):
                if self.h5.isAnalyzed(i):
                    break
                i -= 1
            self.start_frame = i

        # Otherwise start from the beginning.
        else:
            self.h5 = saH5Py.SAH5Py(filename = self.filename,
                                    is_existing = False,
                                    sa_type = sa_type)
            
            # Save analysis parameters.
            etree = parameters.toXMLElementTree(False)
            if (sys.version_info > (3, 0)):
                self.h5.addMetadata(ElementTree.tostring(etree, 'unicode'))
            else:
                self.h5.addMetadata(ElementTree.tostring(etree, 'ISO-8859-1'))

            # Save pixel size.
            self.h5.setPixelSize(parameters.getAttr("pixel_size"))
            
            # Adjust starting frame based on channel 0 offset.
            if (self.offsets[0] != 0):
                assert(self.offsets[0] > 0), "Channel 0 offset cannot be negative."
                self.start_frame = self.offsets[0]
                print("Adjusted start frame to", self.start_frame, "based on channel 0 offset.")

        self.h5.setAnalysisFinished(False)

    def addPeaks(self, peaks, movie_reader):
        assert(len(peaks) == self.n_planes)
        super(MPDataWriter, self).addPeaks(peaks[0], movie_reader)

        if not self.movie_info_set:
            self.h5.addMovieInformation(movie_reader)
            self.movie_info_set = True
        
        for i in range(len(peaks)):
            self.h5.addLocalizations(peaks[i],
                                     movie_reader.getCurrentFrameNumber(),
                                     channel = i)

    def close(self, finished):
        self.h5.setAnalysisFinished(finished)
        self.h5.close(verbose = True)

    
class MPMovieReader(object):
    """
    analysisIO.MovieReader like object for multi-plane data. This is
    primarily designed to be used in the standard storm-analysis 
    analysis pipeline but has some additional functionality to make
    it easier to use in other modules that need to be able to manipulate
    these sorts of movies.
    """
    def __init__(self, base_name = None, parameters = None, **kwds):
        super(MPMovieReader, self).__init__(**kwds)

        self.backgrounds = []
        self.bg_estimators = []
        self.cur_frame = 0
        self.frames = []
        self.max_frame = 0
        self.offsets = []
        self.parameters = parameters
        self.planes = []

        #
        # Load the movies and offsets for each plane/channel. At present
        # multiplane expects the sCMOS camera calibration data.
        #
        calib_name = mpUtil.getCalibrationAttrs(parameters)
        for i, ext in enumerate(mpUtil.getExtAttrs(parameters)):
            movie_name = base_name + parameters.getAttr(ext)
            self.planes.append(analysisIO.FrameReaderSCMOS(parameters = parameters,
                                                           movie_file = movie_name,
                                                           calibration_file = parameters.getAttr(calib_name[i])))

        for offset in mpUtil.getOffsetAttrs(parameters):
            self.offsets.append(parameters.getAttr(offset))

        print("Found data for", len(self.planes), "planes.")

        [self.movie_x, self.movie_y, self.movie_l] = self.planes[0].filmSize()
        self.movie_l -= self.offsets[0]

        # Check if the movies for the other channels (adjusted for their offsets)
        # are shorter than the movie for channel 0.
        #
        for i in range(1, len(self.planes)):
            [px, py, pl] = self.planes[1].filmSize()
            pl -= self.offsets[i]
            if (pl < self.movie_l):
                self.movie_l = pl

        # Assert that all the movies are the same size, at least in x,y.
        for i in range(1, len(self.planes)):
            assert(self.movie_x == self.planes[i].filmSize()[0])
            assert(self.movie_y == self.planes[i].filmSize()[1])

    def close(self):
        for plane in self.planes:
            plane.close()
            
    def getBackground(self, plane):
        if (len(self.backgrounds) > 0):
            return self.backgrounds[plane]
        else:
            return None

    def getCurrentFrameNumber(self):
        return self.cur_frame

    def getFilmSize(self):
        return [self.movie_x, self.movie_y, self.movie_l]
    
    def getFrame(self, plane):
        """
        This returns a particular one of the currently loaded frames.
        """
        return self.frames[plane]

    def getFrames(self, frame_number):
        """
        This loads all the frames for the specified frame_number, corrects
        them for gain and offset (self.planes is a list of analysisIO.FrameReader
        objects) and returns them as a list.
        """
        frames = []
        for i, plane in enumerate(self.planes):
            frames.append(plane.loadAFrame(frame_number + self.offsets[i]))
        return frames

    def getMovieL(self):
        return self.movie_l
    
    def getMovieX(self):
        return self.movie_x

    def getMovieY(self):
        return self.movie_y

    def getNPlanes(self):
        return len(self.planes)
    
    def hashID(self):
        return self.planes[0].hashID()

    def nextFrame(self):
        self.cur_frame += 1
        if (self.cur_frame < self.max_frame):

            # Update background estimate.
            self.backgrounds = []
            for i, bg_estimator in enumerate(self.bg_estimators):
                self.backgrounds.append(bg_estimator.estimateBG(self.cur_frame + self.offsets[i]))

            # Load planes & remove all values less than 1.0 as we are doing MLE fitting.
            self.frames = []
            frames = self.getFrames(self.cur_frame)
            for frame in frames:
                mask = (frame < 1.0)
                if numpy.count_nonzero(mask):
                    frame[mask] = 1.0
                self.frames.append(frame)

            return True
        else:
            return False

    def setup(self, start_frame):

        # Figure out where to start.
        self.cur_frame = start_frame
        if self.parameters.hasAttr("start_frame"):
            if (self.parameters.getAttr("start_frame") > self.cur_frame):
                if (self.parameters.getAttr("start_frame") < self.movie_l):
                    self.cur_frame = self.parameters.getAttr("start_frame") - 1
        
        # Figure out where to stop.
        self.max_frame = self.movie_l

        # If the user specified a max frame then just use it and
        # assume that they knew what they were doing.
        if self.parameters.hasAttr("max_frame"):
            if (self.parameters.getAttr("max_frame") > 0):
                if (self.parameters.getAttr("max_frame") < self.movie_l):
                    self.max_frame = self.parameters.getAttr("max_frame")

        # Configure background estimator, if any.
        #
        # FIXME: Use of a background estimator has not been tested.
        #
        if (self.parameters.getAttr("static_background_estimate", 0) > 0):
            print("Using static background estimator.")
            s_size = self.parameters.getAttr("static_background_estimate")
            for i in range(len(self.planes)):
                bg_est = static_background.StaticBGEstimator(self.planes[i],
                                                             start_frame = self.cur_frame + self.offsets[i],
                                                             sample_size = s_size)
                self.bg_estimators.append(bg_est)
