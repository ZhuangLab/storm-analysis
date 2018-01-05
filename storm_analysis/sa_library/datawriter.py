#!/usr/bin/env python
"""
Writes dax files or tifffile files. This is mostly used
by the simulator.

Hazen 1/18
"""

import numpy
import os
import tifffile


def inferWriter(filename, width = None, height = None):
    """
    Given a file name this will try to return the appropriate
    writer based on the file extension.
    """
    ext = os.path.splitext(filename)[1]
    if (ext == ".dax"):
        return DaxWriter(filename, width = width, height = height)
    elif (ext == ".tif") or (ext == ".tiff"):
        return TiffWriter(filename, width = width, height = height)
    else:
        print(ext, "is not a recognized file type")
        raise IOError("only .dax and .tif are supported (case sensitive..)")

def dummyDaxFile(name, x_size, y_size):
    ddax = DaxWriter(name, width = x_size, height = y_size)
    frame = numpy.ones((x_size, y_size))
    ddax.addFrame(frame)
    ddax.close()

def singleFrameDax(name, frame):
    [fx, fy] = frame.shape
    dax_file = DaxWriter(name, width = fy, height = fx)
    dax_file.addFrame(frame)
    dax_file.close()


class DaxWriter(object):

    def __init__(self, name, width = None, height = None, **kwds):
        super(DaxWriter, self).__init__(**kwds)
        
        self.name = name
        if len(os.path.dirname(name)) > 0:
            self.root_name = os.path.dirname(name) + "/" + os.path.splitext(os.path.basename(name))[0]
        else:
            self.root_name = os.path.splitext(os.path.basename(name))[0]
        self.w = width
        self.h = height
        self.fp = open(self.name, "wb")
        self.l = 0

    def addFrame(self, frame):
        frame = frame.copy()
        frame[(frame < 0)] = 0
        frame[(frame > 65535)] = 65535
        image16 = frame.astype(numpy.uint16)
        image16.tofile(self.fp)
        self.l += 1
        if (self.w is None) or (self.h is None):
            [self.w, self.h] = frame.shape
        else:
            assert(self.w == frame.shape[0])
            assert(self.h == frame.shape[1])

    def close(self):
        self.fp.close()

        self.w = int(self.w)
        self.h = int(self.h)
        
        inf_fp = open(self.root_name + ".inf", "w")
        inf_fp.write("binning = 1 x 1\n")
        inf_fp.write("data type = 16 bit integers (binary, little endian)\n")
        inf_fp.write("frame dimensions = " + str(self.h) + " x " + str(self.w) + "\n")
        inf_fp.write("number of frames = " + str(self.l) + "\n")
        inf_fp.write("Lock Target = 0.0\n")
        if True:
            inf_fp.write("x_start = 1\n")
            inf_fp.write("x_end = " + str(self.h) + "\n")
            inf_fp.write("y_start = 1\n")
            inf_fp.write("y_end = " + str(self.w) + "\n")
        inf_fp.close()


class TiffWriter(object):

    def __init__(self, filename, width = None, height = None, **kwds):
        super(TiffWriter, self).__init__(**kwds)
        self.h = height
        self.w = width
        self.tif_fp = tifffile.TiffWriter(filename)

    def addFrame(self, frame):
        frame = frame.copy()
        frame[(frame < 0)] = 0
        frame[(frame > 65535)] = 65535
        self.tif_fp.save(frame.astype(numpy.uint16))

        # Enforce that all the frames are the same size.
        if (self.h is None) or (self.w is None):
            [self.w, self.h] = frame.shape
        else:
            assert(self.w == frame.shape[0])
            assert(self.h == frame.shape[1])

    def close(self):
        self.tif_fp.close()
        
