#!/usr/bin/env python
"""
Reads Insight3 format binary molecule lists.
Returns results as a Python list of numpy arrays.

Hazen 4/09
"""

import numpy
import os
import struct

from xml.etree import ElementTree
    
import storm_analysis.sa_library.i3dtype as i3dtype


class I3BadStatusException(Exception):
    pass


#
# Functions
#
def _getV(fp, format, size):
    try:
        return struct.unpack(format, fp.read(size))[0]
    except struct.error:
        return 0

def checkStatus(filename):
    with open(filename, "rb") as fp:
        [frames, molecules, version, status] = readHeader(fp, False)
    return (status == 6) and (molecules >= 0) and (version == "M425")

def loadI3File(filename, verbose = True):
    return loadI3FileNumpy(filename, verbose = verbose)

def loadI3FileNumpy(filename, verbose = True):
    with open(filename, "rb") as fp:

        # Read header.
        #
        # This also moves the file pointer to the start of the
        # (binary) localization data.
        #
        [frames, molecules, version, status] = readHeader(fp, verbose)

        # Check status.
        #
        # If this is not 6 then this file was not closed
        # properly, so we'll just stop here.
        #
        if (status != 6):
            print(filename, "was not closed properly, possibly corrupted.")
            return None

        # Read in the localizations.
        #
        # This could include garbage at the end if the file
        # was created by Insight3 or has meta-data.
        #
        data = numpy.fromfile(fp, dtype=i3dtype.i3DataType())
                    
        # Return only the valid localization data.
        return data[:][0:molecules]

def loadI3GoodOnly(filename, verbose = True):
    return loadI3NumpyGoodOnly(filename, verbose = verbose)

def loadI3Metadata(filename, verbose = True):
    """
    Read the metadata, if available. This will always be XML, usually
    it is the analysis settings file.

    We need to be a little careful because if the file was created
    by Insight3 it will have localization data by frame after
    the "master" list.
    """
    with open(filename, "rb") as fp:
        file_size = os.fstat(fp.fileno()).st_size

        # Read header, this also moves the file pointer to
        # the start of the (binary) localization data.
        [frames, molecules, version, status] = readHeader(fp, False)
        locs_end = 16 + molecules * recordSize() + 4

        # Move to the end of the localization data.
        fp.seek(locs_end)
    
        # Check if there is any meta-data.
        if ((fp.tell()+5) < file_size):

            # Check for "<?xml" tag.
            fp_loc = fp.tell()
            if (_getV(fp, "5s", 5).decode() == "<?xml"):
                if verbose:
                    print("Found metadata.")

                # Reset file pointer and read text.
                fp.seek(locs_end)
                return ElementTree.parse(fp).getroot()

            else:
                if verbose:
                    print("No meta data.")

    return None

def loadI3NumpyGoodOnly(filename, verbose = True):
    """
    This filters out molecules with poor z fits (daostorm "3D" fit mode).
    """
    data = loadI3FileNumpy(filename, verbose = verbose)
    return i3dtype.maskData(data, (data['c'] != 9))
            
def readHeader(fp, verbose):
    version = _getV(fp, "4s", 4).decode()
    frames = _getV(fp, "i", 4)
    status = _getV(fp, "i", 4)
    molecules = _getV(fp, "i", 4)
    if verbose:
        print("Version:", version)
        print("Frames:", frames)
        print("Status:", status)
        print("Molecules:", molecules)
        print("")
    return [frames, molecules, version, status]

def recordSize():
    return 4 * i3dtype.getI3DataTypeSize()


class I3Reader(object):
    """
    Binary file reader class.
    """
    def __init__(self, filename, max_to_load = 2000000):
        self.cur_molecule = 0
        self.filename = filename
        self.fp = open(filename, "rb")
        self.localizations = False
        self.record_size = recordSize()
        
        # Load header data
        header_data = readHeader(self.fp, True)
        self.frames = header_data[0]
        self.molecules = header_data[1]
        self.version = header_data[2]
        self.status = header_data[3]

        if (self.status != 6):
            raise I3BadStatusException(filename + " was not closed properly, possibly corrupted.")

        # If the file is small enough, just load all the molecules into memory.
        if (self.molecules < max_to_load):
            self.localizations = loadI3FileNumpy(filename, verbose = False)
            assert (self.molecules == self.localizations.size), "The number of localizations in the file does not match the value in the header."

    def __enter__(self):
        return self

    def __exit__(self, etype, value, traceback):
        if self.fp:
            self.fp.close()

    def close(self):
        self.fp.close()
        
    def getFilename(self):
        return self.filename

    def getMolecule(self, molecule):
        if(molecule < self.molecules):
            if isinstance(self.localizations, numpy.ndarray):
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
        if isinstance(self.localizations, numpy.ndarray):
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
        #
        # FIXME: This is going to return the wrong answer unless
        #        all of the localizations were loaded.
        #
        if(self.molecules > 0):
            mol = self.getMolecule(self.molecules-1)
            return int(mol['fr'][0])
        else:
            return 0

    def getNumberMolecules(self):
        return self.molecules
        
    def findFrame(self, frame):
        def getFrame(molecule):
            mol = self.getMolecule(molecule)
            return int(mol['fr'][0])
        def bin_search(low, high):
            if((low+1)==high):
                return high
            cur = int((low+high)/2)
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

        # Check if we have read all the molecules.
        if(self.cur_molecule>=self.molecules):
            return False

        size = block_size
        
        # Adjust size if we'll read past the end of the file.
        if((self.cur_molecule+size)>self.molecules):
            size = self.molecules - self.cur_molecule

        self.cur_molecule += size

        # Read the data.
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
if (__name__ == "__main__"):
    
    import sys

    from xml.dom import minidom
    from xml.etree import ElementTree

    if (len(sys.argv) != 2):
        print("usage: <bin_file>")
        exit()
    
    print("Checking for meta data.")
    metadata = loadI3Metadata(sys.argv[1], verbose = True)
    if metadata is not None:

        print(" meta data:")
        print("  movie:")
        for node in sorted(metadata.find("movie"), key = lambda node: node.tag):
            print("    " + node.tag.strip() + " - " + node.text.strip())
        print("")
        print("  settings:")
        for node in sorted(metadata.find("settings"), key = lambda node: node.tag):
            print("    " + node.tag.strip() + " - " + node.text.strip())
        print()
    print()

    print("Localization statistics")
    with I3Reader(sys.argv[1]) as i3_in:

        print("Frames:", i3_in.getNumberFrames())

        data = i3_in.nextBlock(good_only = False)
        if data is not False:
            for field in data.dtype.names:
                print(" ", field,"\t",  numpy.mean(data[field]), numpy.std(data[field]), numpy.min(data[field]), numpy.max(data[field]))

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
