#!/usr/bin/python
#
# Read CS high res format file.
#
# Hazen 10/12
#

import numpy
import struct


#
# Functions
#
def _getV(fp, format, size):
    return struct.unpack(format, fp.read(size))[0]

def hrDataType():
    return numpy.dtype([('fr', numpy.int32),    # frame number
                        ('i', numpy.int32),     # object offset
                        ('z', numpy.float32)])  # object intensity

def loadHRFile(filename, verbose = 1):
    fp = open(filename, "rb")

    # Read header
    [x_size, y_size] = readHeader(fp, verbose)

    # Read data
    data = numpy.fromfile(fp, dtype=hrDataType())
    fp.close()

    return [x_size, y_size, data]

def readHeader(fp, verbose):
    x_size = _getV(fp, "i", 4)
    y_size = _getV(fp, "i", 4)
    if verbose:
        print "Size:", x_size, "x", y_size

    # move file pointer to start of the data.
    fp.seek(100)

    return [x_size, y_size]


#
# High res file reader class
#
class HResFile:

    def __init__(self, filename):
        [self.x_size, self.y_size, self.data] = loadHRFile(filename)
        self.first_frame = numpy.min(self.data['fr'])
        self.last_frame = numpy.max(self.data['fr'])

    def getFrame(self, frame_number, binning = 1):
        if (frame_number >= self.first_frame) and (frame_number <= self.last_frame):
            mask = (self.data['fr'] == frame_number)
            
            offset = self.data['i'][mask]
            data = self.data['z'][mask]

            y_vals = numpy.rint(offset/self.x_size).astype(numpy.int32)/binning
            x_vals = numpy.rint(offset%self.x_size).astype(numpy.int32)/binning

            image = numpy.zeros((self.y_size/binning, self.x_size/binning))
            
            for i in range(x_vals.size):
                image[y_vals[i],x_vals[i]] += data[i]

            return image
        else:
            print "Frame number", frame_number, "is out of range (", self.first_frame, ",", self.last_frame, ")"
            return False

    def getSize(self):
        return [self.x_size, self.y_size, self.first_frame, self.last_frame]

    # frame range is inclusive.
    def sumFrames(self, first_frame = -1, last_frame = -1, binning = 1, verbose = False):
        if (first_frame < 0):
            first_frame = self.first_frame
        if (last_frame < 0):
            last_frame = self.last_frame

        image = numpy.zeros((self.y_size,self.x_size))
        for i in range(first_frame, last_frame + 1):
            if verbose and ((i%10) == 0):
                print "Loading frame", i
            image += self.getFrame(i, binning)

        return image
            

if __name__ == "__main__":

    import sys

    # Create a dax movie from a hres file.
    if 0:
        import sa_library.daxwriter as daxwriter

        if (len(sys.argv) != 4):
            print "usage: <in_hres> <out_dax> <binning>"
            exit()

        print "Loading High Res Data"
        hresf = HResFile(sys.argv[1])

        print "Creating Dax File"
        print "  Size info:", hresf.getSize()
        [xs, ys, ff, lf] = hresf.getSize()

        dax_data = daxwriter.DaxWriter(sys.argv[2],0,0)
        binning = int(sys.argv[3])

        for i in range(ff,lf+1):
            print "Creating frame:",i
            frame = hresf.getFrame(i, binning)
            dax_data.addFrame(frame)

        dax_data.close()

    # Create an image from a hres file
    if 1:
        import os

        import sa_library.arraytoimage as arraytoimage
        import sa_library.daxwriter as daxwriter

        if (len(sys.argv) != 3):
            print "usage: <in_hres> <out_img>"
            exit()

        hres = HResFile(sys.argv[1])
        image = hres.sumFrames(verbose = True)

        ext = os.path.splitext(sys.argv[2])[1]
        if (ext == ".png"):
            arraytoimage.singleColorImage(image, sys.argv[2])
        elif (ext == ".dax"):
            daxwriter.singleFrameDax(sys.argv[2], image)
        else:
            print "unrecognized extension ", ext

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
