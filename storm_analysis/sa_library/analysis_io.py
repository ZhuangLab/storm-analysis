#/usr/bin/env python
"""
Input/Output classes that are used by 3D-DAOSTORM, sCMOS, Spliner
and Multiplane analysis.

Hazen 09/17
"""
import numpy
import os
import sys

from xml.etree import ElementTree

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_library.static_background as static_background
import storm_analysis.sa_library.writeinsight3 as writeinsight3


class AnalysisIOException(Exception):
    pass


def isLatestFormat(filename):
    """
    Returns True if the calibration file is the latest format. This
    is mostly used for testing.
    """
    data = numpy.load(filename, allow_pickle = True)
    return (len(data) == 5) and (data[4] == 2)


def loadCMOSCalibration(filename, verbose = False):
    """
    CMOS calibration file reader. This will return the CMOS calibration 
    data as a list [offset, variance, gain, rqe].

    offset - Pixel offset in units of ADU.
    variance - Pixel variance in units of ADU*ADU.
    gain - Pixel gain in units of ADU / photo-electron (e-).
    rqe - Pixel relative quantum efficiency, dimensionless, it should be around 1.0.
   
    filename - The name of calibration file. This should have been saved
               using numpy.save.
    """
    data = numpy.load(filename, allow_pickle = True)
    
    # Check for v0 format.    
    if (len(data) == 3):
        if verbose:
            print("Version 0 sCMOS calibration file detected! Transposing!")

        rqe = numpy.ones_like(data[0])
        data = list(data) + [rqe]
        return map(numpy.transpose, data)

    # Check for v1 format.
    elif (len(data) == 4):
        if verbose:
            print("Version 1 sCMOS calibration file detected! Using 1.0 for relative QE array!")
            
        # v1 format.
        if (data[3] == 1):
            rqe = numpy.ones_like(data[0])
            temp = list(data[:3])
            temp.append(rqe)
            return temp

    # Check for v2 format.
    elif (len(data) == 5):
        
        # v2 format.
        if (data[4] == 2):
            return list(data[:4])

    raise AnalysisIOException("Unknown sCMOS data format.")

    
class DataWriter(object):
    """
    Encapsulate saving the output of the peak finder/fitter.
    """
    def __init__(self, data_file = None, **kwds):
        super(DataWriter, self).__init__(**kwds)

        self.filename = data_file
        self.n_added = 0
        self.start_frame = -1
        self.total_peaks = 0

    def addPeaks(self, peaks, movie_reader):
        self.n_added = peaks["x"].size
        self.total_peaks += self.n_added

    def getNumberAdded(self):
        return self.n_added
    
    def getFilename(self):
        return self.filename

    def getStartFrame(self):
        return self.start_frame
        
    def getTotalPeaks(self):
        return self.total_peaks


class DataWriterHDF5(DataWriter):
    """
    Encapsulate saving the output of the peak finder/fitter to a HDF5 file.
    """
    def __init__(self, parameters = None, sa_type = None, **kwds):
        super(DataWriterHDF5, self).__init__(**kwds)

        self.movie_info_set = False
                
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

        self.h5.setAnalysisFinished(False)

    def addPeaks(self, peaks, movie_reader):
        super(DataWriterHDF5, self).addPeaks(peaks, movie_reader)

        if not self.movie_info_set:
            self.h5.addMovieInformation(movie_reader)
            self.movie_info_set = True

        self.h5.addLocalizations(peaks, movie_reader.getCurrentFrameNumber())

    def close(self, finished):
        self.h5.setAnalysisFinished(finished)
        self.h5.close(verbose = True)

        
class DataWriterI3(DataWriter):
    """
    Encapsulate saving the output of the peak finder/fitter to an Insight3 format file.

    Note: This format is deprecated.
    """
    def __init__(self, parameters = None, **kwds):
        super(DataWriterI3, self).__init__(**kwds)

        raise Exception("Using the Insight3 format for analysis is deprecated!")

        self.pixel_size = parameters.getAttr("pixel_size")
                
        #
        # If the i3 file already exists, read it in, write it
        # out to prepare for starting the analysis from the
        # end of what currently exists.
        #
        # FIXME: If the existing file is really large there
        #        could be problems here as we're going to load
        #        the whole thing into memory.
        #
        if(os.path.exists(self.filename)):
            print("Found", self.filename)
            i3data_in = readinsight3.loadI3File(self.filename)
            if (i3data_in is None) or (i3data_in.size == 0):
                self.start_frame = 0
            else:
                self.start_frame = int(numpy.max(i3data_in['fr']))

            print(" Starting analysis at frame:", self.start_frame)
            self.i3data = writeinsight3.I3Writer(self.filename)
            if (self.start_frame > 0):
                self.i3data.addMolecules(i3data_in)
                self.total_peaks = i3data_in['x'].size
        else:
            self.start_frame = 0
            self.i3data = writeinsight3.I3Writer(self.filename)

    def addPeaks(self, peaks, movie_reader):
        super(DataWriterI3, self).addPeaks(peaks, movie_reader)
        self.i3data.addMultiFitMolecules(peaks,
                                         movie_reader.getMovieX(),
                                         movie_reader.getMovieY(),
                                         movie_reader.getCurrentFrameNumber(),
                                         self.pixel_size)

    def close(self, metadata = None):
        if metadata is None:
            self.i3data.close()
        else:
            self.i3data.closeWithMetadata(metadata)


class FrameReader(object):
    """
    Wraps datareader.Reader, converts frames from ADU to photo-electrons.
    """
    def __init__(self, movie_file = None, parameters = None, **kwds):
        super(FrameReader, self).__init__(**kwds)

        self.gain = None
        self.offset = None
        self.parameters = parameters
        self.rqe = 1.0
        self.verbose = 1
        if self.parameters is not None:
            self.verbose = (self.parameters.getAttr("verbosity") == 1)
        self.movie_data = datareader.inferReader(movie_file)

    def close(self):
        self.movie_data.close()

    def filmSize(self):
        return self.movie_data.filmSize()

    def hashID(self):
        return self.movie_data.hashID()

    def loadAFrame(self, frame_number):

        # Load frame.
        frame = self.movie_data.loadAFrame(frame_number)

        # Convert from ADU to photo-electrons and correct for RQE.
        frame = (frame - self.offset) * (self.gain * self.rqe)

        # Set all values less than 1.0 to 1.0 as we are doing MLE fitting which
        # has zero tolerance for negative numbers..
        #
        mask = (frame < 1.0)
        if (numpy.sum(mask) > 0):
            if self.verbose:
                print(" Removing values < 1.0 in frame {0:0d}".format(frame_number))
            frame[mask] = 1.0

        return frame


class FrameReaderStd(FrameReader):
    """
    Read frames from a 'standard' (as opposed to sCMOS) camera.

    Note: Gain is in units of ADU / photo-electrons.
    """
    def __init__(self, camera_gain = None, camera_offset = None, **kwds):
        super(FrameReaderStd, self).__init__(**kwds)

        if camera_gain is None:
            self.gain = 1.0/self.parameters.getAttr("camera_gain")
            self.offset = self.parameters.getAttr("camera_offset")
        else:
            self.gain = 1.0/camera_gain
            self.offset = camera_offset

        
class FrameReaderSCMOS(FrameReader):
    """
    Read frames from a sCMOS camera.

    Note: Gain is in units of ADU / photo-electrons.
    """
    def __init__(self, calibration_file = None, **kwds):
        super(FrameReaderSCMOS, self).__init__(**kwds)
        
        if calibration_file is None:
            [self.offset, variance, gain, rqe] = loadCMOSCalibration(self.parameters.getAttr("camera_calibration"),
                                                                     verbose = True)
        else:
            [self.offset, variance, gain, rqe] = loadCMOSCalibration(calibration_file,
                                                                     verbose = True)
        self.gain = 1.0/gain
        self.rqe = 1.0/rqe

    
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
        self.cur_frame = -1
        self.frame_reader = frame_reader
        self.frame = None
        self.max_frame = None
        [self.movie_x, self.movie_y, self.movie_l] = frame_reader.filmSize()
        self.parameters = parameters

    def close(self):
        self.frame_reader.close()
        
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
        self.cur_frame += 1
        if (self.cur_frame < self.max_frame):

            # Update background estimate.
            if self.bg_estimator is not None:
                self.background = self.bg_estimator.estimateBG(self.cur_frame)

            # Load frame.
            self.frame = self.frame_reader.loadAFrame(self.cur_frame)

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
