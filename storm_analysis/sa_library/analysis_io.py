#/usr/bin/env python
"""
Input/Output classes that are used by 3D-DAOSTORM, sCMOS, Spliner
and Multiplane analysis.

Hazen 09/17
"""
import numpy
import os

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.static_background as static_background
import storm_analysis.sa_library.writeinsight3 as writeinsight3


class DataWriter(object):
    """
    Encapsulate saving the output of the peak finder/fitter.
    """
    def __init__(self, data_file = None, parameters = None, **kwds):
        super(DataWriter, self).__init__(**kwds)

        self.filename = data_file
        self.inverted = (parameters.getAttr("orientation", "normal") == "inverted")
        self.pixel_size = parameters.getAttr("pixel_size")
        self.n_added = 0
                
        #
        # If the i3 file already exists, read it in, write it
        # out to prepare for starting the analysis from the
        # end of what currently exists.
        #
        # FIXME: If the existing file is really large there
        #        could be problems here as we're going to load
        #        the whole thing into memory.
        #
        self.total_peaks = 0
        if(os.path.exists(data_file)):
            print("Found", data_file)
            i3data_in = readinsight3.loadI3File(data_file)
            if i3data_in is None:
                self.start_frame = 0
            else:
                self.start_frame = int(numpy.max(i3data_in['fr']))

            print(" Starting analysis at frame:", self.start_frame)
            self.i3data = writeinsight3.I3Writer(data_file)
            if (self.start_frame > 0):
                self.i3data.addMolecules(i3data_in)
                self.total_peaks = i3data_in['x'].size
        else:
            self.start_frame = 0
            self.i3data = writeinsight3.I3Writer(data_file)

    def addPeaks(self, peaks, movie_reader):
        self.n_added = peaks.shape[0]
        self.i3data.addMultiFitMolecules(peaks,
                                         movie_reader.getMovieX(),
                                         movie_reader.getMovieY(),
                                         movie_reader.getCurrentFrameNumber(),
                                         self.pixel_size,
                                         self.inverted)
        self.total_peaks += peaks.shape[0]

    def close(self, metadata = None):
        if metadata is None:
            self.i3data.close()
        else:
            self.i3data.closeWithMetadata(metadata)

    def getNumberAdded(self):
        return self.n_added
    
    def getFilename(self):
        return self.filename

    def getStartFrame(self):
        return self.start_frame
        
    def getTotalPeaks(self):
        return self.total_peaks
    

class FrameReader(object):
    """
    Wraps datareader.Reader, converts frames from ADU to photo-electrons.
    """
    def __init__(self, movie_file = None, **kwds):
        super(FrameReader, self).__init__(**kwds)
        
        self.movie_data = datareader.inferReader(movie_file)

    def filmSize(self):
        return self.movie_data.filmSize()

    def hashID(self):
        return self.movie_data.hashID()

    def loadAFrame(self, frame_number):

        # Load frame.
        frame = self.movie_data.loadAFrame(frame_number)

        # Convert from ADU to photo-electrons.
        frame = (frame - self.offset) * self.gain

        return frame
        
        
class FrameReaderStd(FrameReader):
    """
    Read frames from a 'standard' (as opposed to sCMOS) camera.

    Note: Gain is in units of ADU / photo-electrons.
    """
    def __init__(self, parameters = None, **kwds):
        super(FrameReaderStd, self).__init__(**kwds)

        self.gain = 1.0/parameters.getAttr("camera_gain")
        self.offset = parameters.getAttr("camera_offset")

        
class FrameReaderSCMOS(FrameReader):
    """
    Read frames from a sCMOS camera.

    Note: Gain is in units of ADU / photo-electrons.
    """
    def __init__(self, parameters = None, **kwds):
        super(FrameReaderSCMOS, self).__init__(**kwds)

        [self.offset, variance, gain] = numpy.load(parameters.getAttr("camera_calibration"))
        self.gain = 1.0/gain

    
class MovieReader(object):
    """
    Encapsulate getting the frames of a movie from a datareader.Reader, handling
    static_background, which frame to start on, etc...

    The primary reason for using a class like this is to add the flexibility
    necessary in order to make the peakFinding() function also work with
    multi-channel data / analysis.
    """
    def __init__(self, frame_reader = None, parameters = None, **kwds):
        super(MovieReader, self).__init__(**kwds)

        self.background = None
        self.bg_estimator = None
        self.cur_frame = 0
        self.frame_reader = frame_reader
        self.frame = None
        self.max_frame = None
        [self.movie_x, self.movie_y, self.movie_l] = frame_reader.filmSize()
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
        return self.frame_reader.hashID()

    def nextFrame(self):
        if (self.cur_frame < self.max_frame):

            # Update background estimate.
            if self.bg_estimator is not None:
                self.background = self.bg_estimator.estimateBG(self.cur_frame)

            # Load frame & remove all values less than 1.0 as we are doing MLE fitting.
            self.frame = self.frame_reader.loadAFrame(self.cur_frame)
            mask = (self.frame < 1.0)
            if (numpy.sum(mask) > 0):
                print(" Removing values < 1.0 in frame", self.cur_frame)
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
            self.bg_estimator = static_background.StaticBGEstimator(self.frame_reader,
                                                                    start_frame = self.cur_frame,
                                                                    sample_size = s_size)
