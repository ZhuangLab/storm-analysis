#!/usr/bin/python
#
# Classes that handles reading STORM movie files. Currently this
# is limited to the .dax and .spe formats
#
# Hazen 05/12
#

import numpy
import os
import re


#
# The superclass containing those functions that 
# are common to reading a STORM movie file.
#
# Subclasses should implement:
#  1. __init__(self, filename, verbose = False)
#     This function should open the file and extract the
#     various key bits of meta-data such as the size in XY
#     and the length of the movie.
#
#  2. loadAFrame(self, frame_number)
#     Load the requested frame and return it as numpy array.
#
class Reader:
    # close the file on cleanup
    def __del__(self):
        if self.fileptr:
            self.fileptr.close()

    # average multiple frames in a movie
    def averageFrames(self, start = False, end = False):
        if (not start):
            start = 0
        if (not end):
            end = self.number_frames 

        length = end - start
        average = numpy.zeros((self.image_width, self.image_height), numpy.float)
        for i in range(length):
            average += self.loadAFrame(i + start)
            
        average = average/float(length)
        return average

    # returns the film name
    def filmFilename(self):
        return self.filename

    # returns the film size
    def filmSize(self):
        return [self.image_width, self.image_height, self.number_frames]

    # returns the picture x,y location, if available
    def filmLocation(self):
        if hasattr(self, "stage_x"):
            return [self.stage_x, self.stage_y]
        else:
            return [0.0, 0.0]

    # returns the film focus lock target
    def lockTarget(self):
        if hasattr(self, "lock_target"):
            return self.lock_target
        else:
            return 0.0

    # returns the scale used to display the film when
    # the picture was taken.
    def filmScale(self):
        if hasattr(self, "scalemin") and hasattr(self, "scalemax"):
            return [self.scalemin, self.scalemax]
        else:
            return [100, 2000]


#
# Dax reader class. This is a Zhuang lab custom format.
#
class DaxReader(Reader):
    # dax specific initialization
    def __init__(self, filename, verbose = 0):
        # save the filenames
        self.filename = filename
        dirname = os.path.dirname(filename)
        if (len(dirname) > 0):
            dirname = dirname + "/"
        self.inf_filename = dirname + os.path.splitext(os.path.basename(filename))[0] + ".inf"

        # defaults
        self.image_height = None
        self.image_width = None

        # extract the movie information from the associated inf file
        size_re = re.compile(r'frame dimensions = ([\d]+) x ([\d]+)')
        length_re = re.compile(r'number of frames = ([\d]+)')
        endian_re = re.compile(r' (big|little) endian')
        stagex_re = re.compile(r'Stage X = ([\d\.\-]+)')
        stagey_re = re.compile(r'Stage Y = ([\d\.\-]+)')
        lock_target_re = re.compile(r'Lock Target = ([\d\.\-]+)')
        scalemax_re = re.compile(r'scalemax = ([\d\.\-]+)')
        scalemin_re = re.compile(r'scalemin = ([\d\.\-]+)')

        inf_file = open(self.inf_filename, "r")
        while 1:
            line = inf_file.readline()
            if not line: break
            m = size_re.match(line)
            if m:
                self.image_height = int(m.group(1))
                self.image_width = int(m.group(2))
            m = length_re.match(line)
            if m:
                self.number_frames = int(m.group(1))
            m = endian_re.search(line)
            if m:
                if m.group(1) == "big":
                    self.bigendian = 1
                else:
                    self.bigendian = 0
            m = stagex_re.match(line)
            if m:
                self.stage_x = float(m.group(1))
            m = stagey_re.match(line)
            if m:
                self.stage_y = float(m.group(1))
            m = lock_target_re.match(line)
            if m:
                self.lock_target = float(m.group(1))
            m = scalemax_re.match(line)
            if m:
                self.scalemax = int(m.group(1))
            m = scalemin_re.match(line)
            if m:
                self.scalemin = int(m.group(1))

        inf_file.close()

        # set defaults, probably correct, but warn the user 
        # that they couldn't be determined from the inf file.
        if not self.image_height:
            print "Could not determine image size, assuming 256x256."
            self.image_height = 256
            self.image_width = 256

        # open the dax file
        if os.path.exists(filename):
            self.fileptr = open(filename, "rb")
        else:
            self.fileptr = 0
            if verbose:
                print "dax data not found", filename

    # load a frame & return it as a numpy array
    def loadAFrame(self, frame_number):
        if self.fileptr:
            assert frame_number >= 0, "frame_number must be greater than or equal to 0"
            assert frame_number < self.number_frames, "frame number must be less than " + str(self.number_frames)
            self.fileptr.seek(frame_number * self.image_height * self.image_width * 2)
            image_data = numpy.fromfile(self.fileptr,dtype='int16', count = self.image_height * self.image_width)
            image_data = numpy.transpose(numpy.reshape(image_data, [self.image_width, self.image_height]))
            if self.bigendian:
                image_data.byteswap(True)
            return image_data

#
# SPE (Roper Scientific) reader class.
#
class SPEReader(Reader):
    # dax specific initialization
    def __init__(self, filename, verbose = 0):
        # save the filename
        self.filename = filename

        # open the file & read the header
        self.header_size = 4100
        self.fileptr = open(filename, "rb")

        self.fileptr.seek(42)
        self.image_width = int(numpy.fromfile(self.fileptr, numpy.uint16, 1)[0])
        self.fileptr.seek(656)
        self.image_height = int(numpy.fromfile(self.fileptr, numpy.uint16, 1)[0])
        self.fileptr.seek(1446)
        self.number_frames = int(numpy.fromfile(self.fileptr, numpy.uint32, 1)[0])

        self.fileptr.seek(108)
        image_mode = int(numpy.fromfile(self.fileptr, numpy.uint16, 1)[0])
        if (image_mode == 0):
            self.image_size = 4 * self.image_width * self.image_height
            self.image_mode = numpy.float32
        elif (image_mode == 1):
            self.image_size = 4 * self.image_width * self.image_height
            self.image_mode = numpy.uint32
        elif (image_mode == 2):
            self.image_size = 2 * self.image_width * self.image_height
            self.image_mode = numpy.int16
        elif (image_mode == 3):
            self.image_size = 2 * self.image_width * self.image_height
            self.image_mode = numpy.uint16
        else:
            print "unrecognized spe image format: ", image_mode

    # load a frame & return it as a numpy array
    def loadAFrame(self, frame_number):
        if self.fileptr:
            assert frame_number >= 0, "frame_number must be greater than or equal to 0"
            assert frame_number < self.number_frames, "frame number must be less than " + str(self.number_frames)
            self.fileptr.seek(self.header_size + frame_number * self.image_size)
            image_data = numpy.fromfile(self.fileptr, dtype=self.image_mode, count = self.image_height * self.image_width)
            image_data = numpy.transpose(numpy.reshape(image_data, [self.image_width, self.image_height]))
            return image_data


#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
