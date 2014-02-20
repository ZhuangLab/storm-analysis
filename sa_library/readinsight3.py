#!/usr/bin/python
#
# Reads Insight3 format binary molecule lists.
# Returns results as a Python list of numpy arrays.
#
# Hazen 4/09
#

import numpy
import struct

import i3dtype


#
# Functions
#
def _getV(fp, format, size):
    return struct.unpack(format, fp.read(size))[0]

def loadI3File(filename, verbose = 1):
    return loadI3FileNumpy(filename, verbose = verbose)

def loadI3FileNumpy(filename, verbose = 1):
    fp = open(filename, "rb")

    # Read header
    [frames, molecules, version, status] = readHeader(fp, verbose)

    data = numpy.fromfile(fp, dtype=i3dtype.i3DataType())

    # If the analysis crashed, the molecule list may still 
    # be valid, but the molecule number will be incorrect.
    if(molecules==0):
        print "File appears empty, trying to load anyway."
        try:
            molecules = numpy.max(data['fr'])
        except:
            molecules = 0
    data = data[:][0:molecules]
    fp.close()
    return data

def loadI3GoodOnly(filename, verbose = 1):
    return loadI3NumpyGoodOnly(filename, verbose = verbose)

def loadI3NumpyGoodOnly(filename, verbose = 1):
    data = loadI3FileNumpy(filename, verbose = verbose)
    return i3dtype.maskData(data, (data['c'] != 9))

def readHeader(fp, verbose):
    version = _getV(fp, "4s", 4)
    frames = _getV(fp, "i", 4)
    status = _getV(fp, "i", 4)
    molecules = _getV(fp, "i", 4)
    if verbose:
        print "Version:", version
        print "Frames:", frames
        print "Status:", status
        print "Molecules:", molecules
        print ""
    return [frames, molecules, version, status]


#
# Bin reader class.
#
class I3Reader:
    def __init__(self, filename, max_to_load = 2000000):
        self.cur_molecule = 0
        self.filename = filename
        self.fp = open(filename, "rb")
        self.localizations = False
        self.record_size = 4 * i3dtype.getI3DataTypeSize()

        # Load header data
        header_data = readHeader(self.fp, 1)
        self.frames = header_data[0]
        self.molecules = header_data[1]
        self.version = header_data[2]
        self.status = header_data[3]

        # If the file is small enough, just load all the molecules into memory.
        if (self.molecules < max_to_load):
            self.localizations = loadI3FileNumpy(filename, verbose = False)

    def close(self):
        self.fp.close()

    def getFilename(self):
        return self.filename

    def getMolecule(self, molecule):
        if(molecule<self.molecules):
            if (type(self.localizations) == type(numpy.array([]))):
                return self.localizations[molecule:molecule+1]
            else:
                cur = self.fp.tell()
                self.fp.seek(16 + molecule*self.record_size)
                data = numpy.fromfile(self.fp,
                                      dtype=i3dtype.i3DataType(),
                                      count=1)
                self.fp.seek(cur)
                return data

    def getMoleculesInFrame(self, frame, good_only = True):
        return self.getMoleculesInFrameRange(frame, frame+1, good_only)

    def getMoleculesInFrameRange(self, start, stop, good_only = True):
        start_mol_num = self.findFrame(start)
        stop_mol_num = self.findFrame(stop)
        if (type(self.localizations) == type(numpy.array([]))):
            data = self.localizations[start_mol_num:stop_mol_num]
        else:
            cur = self.fp.tell()
            self.fp.seek(16 + start_mol_num*self.record_size)
            size = stop_mol_num - start_mol_num
            data = numpy.fromfile(self.fp,
                                  dtype=i3dtype.i3DataType(),
                                  count=size)
            self.fp.seek(cur)  
        if good_only:
            return i3dtype.maskData(data, (data['c'] != 9))
        else:
            return data

    def getNumberFrames(self):
        mol = self.getMolecule(self.molecules-1)
        return int(mol['fr'][0])
        
    def findFrame(self, frame):
        def getFrame(molecule):
            mol = self.getMolecule(molecule)
            return int(mol['fr'][0])
        def bin_search(low, high):
            if((low+1)==high):
                return high
            cur = (low+high)/2
            mid = getFrame(cur)
            if(mid>=frame):
                return bin_search(low,cur)
            elif(mid<frame):
                return bin_search(cur,high)
        if(frame <= getFrame(0)):
            return 0
        elif(frame > getFrame(self.molecules-1)):
            return self.molecules
        else:
            return bin_search(0,self.molecules-1)

    def nextBlock(self, block_size = 400000, good_only = True):

        # check if we have read all the molecules.
        if(self.cur_molecule>=self.molecules):
            return False

        size = block_size
        
        # adjust size if we'll read past the end of the file.
        if((self.cur_molecule+size)>self.molecules):
            size = self.molecules - self.cur_molecule

        self.cur_molecule += size

        # read the data
        data = numpy.fromfile(self.fp,
                              dtype=i3dtype.i3DataType(),
                              count=size)

        if good_only:
            return i3dtype.maskData(data, (data['c'] != 9))
        else:
            return data
        
    def resetFp(self):
        self.fp.seek(16)
        self.cur_molecule = 0


#
# Testing
#
if __name__ == "__main__":
    import sys
    i3_in = I3Reader(sys.argv[1])

    if 1:
        print i3_in.getNumberFrames()

        data = i3_in.nextBlock()
        for field in data.dtype.names:
            print " ", field,"\t",  numpy.mean(data[field]), numpy.std(data[field]), numpy.min(data[field]), numpy.max(data[field])

    if 0:
        print i3_in.getNumberFrames()

        mol = i3_in.getMolecule(int(sys.argv[2]))
        for field in mol.dtype.names:
            print " ", field,"\t",  mol[field][0]

    if 0:
        mol_num = i3_in.findFrame(int(sys.argv[2]))
        print mol_num
        print i3_in.getMolecule(mol_num-1)['fr'][0]
        print i3_in.getMolecule(mol_num)['fr'][0]

    if 0:
        mols = i3_in.getMoleculesInFrameRange(int(sys.argv[2]),
                                              int(sys.argv[3]))
        print numpy.min(mols['fr']), numpy.max(mols['fr'])


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
