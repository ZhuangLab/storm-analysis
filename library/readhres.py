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


    #
    # This does not work, but I think that it should..
    #
    def averageFrames(self, first_frame = -1, last_frame = -1):
        if (first_frame < 0):
            first_frame = self.first_frame
        if (last_frame < 0):
            last_frame = self.last_frame

        mask = (self.data['fr'] >= first_frame) & (self.data['fr'] <= last_frame)

        offset = self.data['i'][mask]
        data = self.data['z'][mask]

        y_vals = numpy.rint(offset/self.x_size).astype(numpy.int32)
        x_vals = numpy.rint(offset%self.x_size).astype(numpy.int32)
        
        image = numpy.zeros((self.y_size,self.x_size))
        
        image[y_vals,x_vals] += data

        return image
            
    def getFrame(self, frame_number):
        if (frame_number >= self.first_frame) and (frame_number <= self.last_frame):
            mask = (self.data['fr'] == frame_number)
            
            offset = self.data['i'][mask]
            data = self.data['z'][mask]

            y_vals = numpy.rint(offset/self.x_size).astype(numpy.int32)
            x_vals = numpy.rint(offset%self.x_size).astype(numpy.int32)

            image = numpy.zeros((self.y_size,self.x_size))
            
            image[y_vals,x_vals] = data

            return image
        else:
            print "Frame number", frame_number, "is out of range (", self.first_frame, ",", self.last_frame, ")"
            return False

    def getSize(self):
        return [self.x_size, self.y_size, self.first_frame, self.last_frame]


if __name__ == "__main__":

    import sys

    import library.daxwriter as daxwriter

    print "Loading High Res Data"
    hresf = HResFile(sys.argv[1])

    print "Creating Dax File"
    print "  Size info:", hresf.getSize()
    [xs, ys, ff, lf] = hresf.getSize()

    dax_data = daxwriter.DaxWriter(sys.argv[2],0,0)

    for i in range(ff,lf+1):
        print "Creating frame:",i
        frame = hresf.getFrame(i)
        dax_data.addFrame(frame)

    dax_data.close()

