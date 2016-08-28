#!/usr/bin/python
#
# Writes dax files. Designed for use by analysis & simulation
# not so much for use by data taking code.
#
# Hazen 4/09
#

import numpy
import os

class DaxWriter():

    def __init__(self, name, w, h):
        self.name = name
        if len(os.path.dirname(name)) > 0:
            self.root_name = os.path.dirname(name) + "/" + os.path.splitext(os.path.basename(name))[0]
        else:
            self.root_name = os.path.splitext(os.path.basename(name))[0]
        self.w = int(w)
        self.h = int(h)
        self.fp = open(self.name, "wb")
        self.l = 0

    def addFrame(self, frame):
        frame = frame.copy()
        mask = (frame < 0)
        frame[mask] = 0
        mask = (frame > 65535)
        frame[mask] = 65535
        image16 = frame.astype('Int16')
        image16 = numpy.transpose(image16)
        image16.byteswap(True)
        image16.tofile(self.fp)
        self.l += 1
        if (self.w == 0) or (self.h == 0):
            [self.w, self.h] = frame.shape

    def close(self):
        self.fp.close()
        inf_fp = open(self.root_name + ".inf", "w")
        inf_fp.write("binning = 1 x 1\n")
        inf_fp.write("data type = 16 bit integers (binary, big endian)\n")
        inf_fp.write("frame dimensions = " + str(self.w) + " x " + str(self.h) + "\n")
        inf_fp.write("number of frames = " + str(self.l) + "\n")
        inf_fp.write("Lock Target = 0.0\n")
        if 1:
            inf_fp.write("x_start = 1\n")
            inf_fp.write("x_end = " + str(self.w) + "\n")
            inf_fp.write("y_start = 1\n")
            inf_fp.write("y_end = " + str(self.h) + "\n")
        inf_fp.close()

def dummyDaxFile(name, x_size, y_size):
    ddax = DaxWriter(name, x_size, y_size)
    frame = numpy.ones((x_size, y_size))
    ddax.addFrame(frame)
    ddax.close()

def singleFrameDax(name, frame):
    [fx, fy] = frame.shape
    dax_file = DaxWriter(name, fy, fx)
    dax_file.addFrame(frame)
    dax_file.close()

