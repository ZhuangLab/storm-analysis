#!/usr/bin/env python
"""
Analysis IO specialized for multiplane fitting.

Hazen 09/17
"""
import numpy

import storm_analysis.multi_plane.mp_utilities as mpUtil

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.writeinsight3 as writeinsight3


class MPDataWriter(analysisIO.DataWriter):
    """
    Data writer specialized for multi-plane data.
    """
    def __init__(self, parameters = None, **kwds):
        super(MPDataWriter, self).__init__(**kwds)

        self.pixel_size = parameters.getAttr("pixel_size")
                
        self.offsets = []

        # Figure out how many planes there are.
        self.n_planes = len(mpUtil.getExtAttrs(parameters))

        # Save frame offsets for each plane.
        for offset in mpUtil.getOffsetAttrs(parameters):
            self.offsets.append(parameters.getAttr(offset))        

        # Adjust starting frame based on channel 0 offset.
        if (self.start_frame > 0) and (self.offsets[0] != 0):
            self.start_frame += self.offsets[0]
            print("Adjusted start frame to", self.start_frame, "based on channel 0 offset.")
        
        # Create writers.
        assert(self.start_frame == 0)
        self.i3_writers = [writeinsight3.I3Writer(self.filename)]
        for i in range(1, self.n_planes):
            fname = self.filename[:-4] + "_ch" + str(i) + ".bin"
            self.i3_writers.append(writeinsight3.I3Writer(fname))

    def addPeaks(self, peaks, movie_reader):
        assert(len(peaks) == self.n_planes)

        for i in range(len(peaks)):
            self.i3_writers[i].addMultiFitMolecules(peaks[i],
                                                    movie_reader.getMovieX(),
                                                    movie_reader.getMovieY(),
                                                    movie_reader.getCurrentFrameNumber() + self.offsets[i],
                                                    self.pixel_size)

        self.n_added = peaks[0]["x"].size
        self.total_peaks += self.n_added

    def close(self, metadata = None):
        for i3w in self.i3_writers:
            if metadata is None:
                i3w.close()
            else:
                i3w.closeWithMetadata(metadata)

    
class MPMovieReader(object):
    """
    analysisIO.MovieReader like object for multi-plane data.

    Note: This uses channel 0 as the reference length and assumes
          that the movies for all the other channels are at least
          this length or longer.
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
            self.planes.append(analysisIO.FrameReaderSCMOS(movie_file = movie_name,
                                                           calibration_file = parameters.getAttr(calib_name[i])))

        for offset in mpUtil.getOffsetAttrs(parameters):
            self.offsets.append(parameters.getAttr(offset))

        print("Found data for", len(self.planes), "planes.")

        [self.movie_x, self.movie_y, self.movie_l] = self.planes[0].filmSize()

        # Assert that all the movies are the same size, at least in x,y.
        for i in range(1, len(self.planes)):
            assert(self.movie_x == self.planes[i].filmSize()[0])
            assert(self.movie_y == self.planes[i].filmSize()[1])
        
    def getBackground(self, plane):
        if (len(self.backgrounds) > 0):
            return self.backgrounds[plane]
        else:
            return None

    def getCurrentFrameNumber(self):
        return self.cur_frame
    
    def getFrame(self, plane):
        return self.frames[plane]

    def getMovieL(self):
        return self.movie_l
    
    def getMovieX(self):
        return self.movie_x

    def getMovieY(self):
        return self.movie_y
    
    def hashID(self):
        return self.planes[0].hashID()

    def nextFrame(self):
        if (self.cur_frame < self.max_frame):

            # Update background estimate.
            self.backgrounds = []
            for i, bg_estimator in enumerate(self.bg_estimators):
                self.backgrounds.append(bg_estimator.estimateBG(self.cur_frame + self.offsets[i]))

            # Load planes & remove all values less than 1.0 as we are doing MLE fitting.
            self.frames = []
            for i, plane in enumerate(self.planes):
                frame = plane.loadAFrame(self.cur_frame + self.offsets[i])
                mask = (frame < 1.0)
                if (numpy.sum(mask) > 0):
                    frame[mask] = 1.0
                self.frames.append(frame)

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

        # Adjust movie length based on channel 0 offset, if any.
        #
        # FIXME: Need to also query other channels in case they are shorter.
        #
        if (len(self.offsets) > 0):
            self.movie_l -= self.offsets[0]

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
