#!/usr/bin/python
#
# Populates 2D or 3D arrays with data from an Insight3 analysis file.
#
# Hazen 7/09
#
# Modified (incompletely) to support lazy loading, which helps
# when dealing with super huge insight3 files.
#
# Hazen 11/11
# 

import numpy
import os
import sys

import datareader
import grid_c
import i3dtype
import readinsight3
import regfilereader


#
# Determine the film size.
#
def getFilmSize(filename, i3_data):
    names = [filename[:-9], filename[:-10]]
    extensions = [".dax", ".spe", ".tif"]
    for name in names:
        for ext in extensions:
            if os.path.exists(name + ext):
                movie_file = datareader.inferReader(name + ext)
                return movie_file.filmSize()

    film_l = int(numpy.max(i3_data['fr']))+1
    print "Could not find movie file for", filename, "assuming 256x256x" + str(film_l)
    return [256, 256, film_l]


#
# Generic Insight3 grid class.
#
class I3GGeneric():
    def __init__(self, filename, scale = 4, verbose = True):

        # Setup names.
        self.filename = filename
        self.fullname = filename
        self.dirname = os.path.dirname(filename)
        if (len(self.dirname) > 0):
            self.dirname = self.dirname + "/"
            self.filename = os.path.splitext(os.path.basename(filename))[0]
            
        # Other class variables
        self.scale = int(scale)
        self.z_range = 1000.0


#
# The I3 grid class.
#
# This class will attempt to load the entire localization list
# into memory. This is fine for smaller data sets but can
# be problematic for large data sets.
#
class I3GData(I3GGeneric):
    def __init__(self, filename, scale = 4, verbose = True):
        I3GGeneric.__init__(self,
                            filename,
                            scale = scale,
                            verbose = verbose)
        
        self.i3data = readinsight3.loadI3GoodOnly(filename, verbose = verbose)
        self.i3data['fr'] -= 1

        # Determine film size.
        [image_x, image_y, self.film_l] = getFilmSize(filename, self.i3data)
        self.im_size = [image_x, image_y]

        # Determine what channels the image has.
        self.channels = []
        for i in range(10):
            mask = (self.i3data['c'] == i)
            if mask.sum() > 0:
                self.channels.append(i)

    # Utility
    def applyXYDriftCorrection(self, dx, dy):
        if(type(dx)==type(numpy.array([]))):
            f = self.i3data['fr']
            i = numpy.arange(f.size, dtype=int)
            dx = dx.astype(numpy.float32)
            dy = dy.astype(numpy.float32)
            self.i3data['xc'][i] = self.i3data['x'] + dx[f]
            self.i3data['yc'][i] = self.i3data['y'] + dy[f]
        else:
            self.i3data['xc'] = self.i3data['x'] + dx
            self.i3data['yc'] = self.i3data['y'] + dy

    def applyZDriftCorrection(self, dz):
        if(type(dz)==type(numpy.array([]))):
            f = self.i3data['fr']
            i = numpy.arange(f.size, dtype=int)
            self.i3data['zc'][i] = self.i3data['z'] + dz[f]
        else:
            self.i3data['zc'] = self.i3data['z'] + dz

    def get2D(self, zmin = -1000.0, zmax = 1000.0):
        return self.i3To2DGridAllChannelsMerged(verbose = 0, zmin = zmin, zmax = zmax)

    def get3D(self, zbins, zmin = -1000.0, zmax = 1000.0):
        return self.i3To3DGridAllChannelsMerged(zbins, zmin = zmin, zmax = zmax, verbose = 0)

    def getData(self):
        return self.data

    def getDirname(self):
        return self.dirname

    def getFilename(self):
        return self.filename

    def getFilmLength(self):
        return self.film_l

    def getFullname(self):
        return self.fullname

    def getImageSize(self):
        return self.im_size

    def getNumberMolecules(self):
        return self.i3data['x'].size

    def getScale(self):
        return self.scale

    def getXY(self):
        return [self.i3data['xc'],
                self.i3data['yc']]

    def getXYZ(self):
        return [self.i3data['xc'],
                self.i3data['yc'],
                self.i3data['zc']]

    def getXYZCat(self):
        return [self.i3data['xc'],
                self.i3data['yc'],
                self.i3data['zc'],
                self.i3data['c']]

    def getXYZICat(self):
        return [self.i3data['xc'],
                self.i3data['yc'],
                self.i3data['zc'],
                self.i3data['i'],
                self.i3data['c']]

    def getZRange(self):
        return self.z_range

    def offsetX(self, dx):
        self.i3data['xc'] += dx

    def offsetY(self, dy):
        self.i3data['yc'] += dy

    def offsetZ(self, dz):
        self.i3data['zc'] += dz

    def setScale(self, scale):
        self.scale = scale

    # Gridding data
    def i3To2DGrid(self, fmin = 0, fmax = 500000, zmin = -1000.0, zmax = 1000.0, uncorrected = False, matrix = False, translate = False, verbose = True):

        if uncorrected:
            [x, y, z] = [self.i3data['x'],
                         self.i3data['y'],
                         self.i3data['z']]
        else:
            [x, y, z] = [self.i3data['xc'],
                         self.i3data['yc'],
                         self.i3data['zc']]
            
        cat = self.i3data['c']
        f = self.i3data['fr']

        [image_x, image_y] = self.im_size
        scale = int(self.scale)

        if type(matrix) == type(numpy.array([])):
            if translate:
                [x, y] = regfilereader.applyTransform(matrix, x, y)
            else:
                [x, y] = regfilereader.applyTransformNoTranslation(matrix, x, y)

        max_max = 0.0
        max_counts = []
        image_data = []
        for channel in self.channels:
            mask = (cat == channel) & (f >= fmin) & (f < fmax) & (z > zmin) & (z < zmax)
            #print numpy.sum(mask)
            i_x = numpy.round(x[mask] * scale).astype(int)
            i_y = numpy.round(y[mask] * scale).astype(int)
            if (i_x.shape[0] > 0):
                channel_data = grid_c.grid2D(i_x,i_y,(image_x*scale,image_y*scale))
            else:
                channel_data = numpy.zeros((image_x*scale, image_y*scale))
            image_data.append(channel_data.astype(numpy.float32))
            max_count = numpy.max(channel_data)
            if max_count > max_max:
                max_max = max_count
            max_counts.append(max_count)

        return [image_data, self.channels, max_counts, max_max]

    def i3To2DGridAllChannelsMerged(self, fmin = 0, fmax = 500000, zmin = -1000.0, zmax = 1000.0, uncorrected = False, verbose = True):
        [image_data, channels, max_counts, max_max] = self.i3To2DGrid(fmin = fmin, 
                                                                      zmin = zmin,
                                                                      zmax = zmax,
                                                                      fmax = fmax, 
                                                                      uncorrected = uncorrected,
                                                                      verbose = verbose)
        if (len(image_data) > 0):
            merged_image = image_data[0]
            for i in range(len(image_data)-1):
                merged_image += image_data[i+1]
            return merged_image
        else:
            return numpy.ones(100,100)

    def i3To3DGrid(self, z_bins, fmin = 0, fmax = 500000, zmin = -1000.0, zmax = 1000.0, uncorrected = False, verbose = True):

        if uncorrected:
            [x, y, z] = [self.i3data['x'],
                         self.i3data['y'],
                         self.i3data['z']]
        else:
            [x, y, z] = [self.i3data['xc'],
                         self.i3data['yc'],
                         self.i3data['zc']]
            
        cat = self.i3data['c']
        f = self.i3data['fr']

        [image_x, image_y] = self.im_size
        xy_scale = int(self.scale)
        z_bins = int(z_bins)

        max_max = 0.0
        max_counts = []
        image_data = []
        for channel in self.channels:
            mask = (cat == channel) & (f >= fmin) & (f < fmax) & (z > zmin) & (z < zmax)
            i_x = numpy.round(x[mask] * xy_scale).astype(int)
            i_y = numpy.round(y[mask] * xy_scale).astype(int)
            i_z = numpy.round((z[mask] + 500.0) * float(z_bins)/1000.0).astype(int)
            if (i_x.shape[0] > 0):
                channel_data = grid_c.grid3D(i_x, i_y, i_z, (image_x*xy_scale, image_y*xy_scale, z_bins))
            else:
                channel_data = numpy.zeros((image_x*xy_scale, image_y*xy_scale, z_bins))
            image_data.append(channel_data.astype(numpy.float32))
            max_count = numpy.max(channel_data)
            if max_count > max_max:
                max_max = max_count
            max_counts.append(max_count)

        return [image_data, self.channels, max_counts, max_max]

    def i3To3DGridAllChannelsMerged(self, z_bins, fmin = 0, fmax = 500000, zmin = -1000.0, zmax = 1000.0, uncorrected = False, verbose = True):
        [image_data, channels, max_counts, max_max] = self.i3To3DGrid(z_bins, 
                                                                      fmin = fmin, 
                                                                      fmax = fmax,
                                                                      zmin = zmin,
                                                                      zmax = zmax,
                                                                      uncorrected = uncorrected,
                                                                      verbose = verbose)
        merged_image = image_data[0]
        for i in range(len(image_data)-1):
            merged_image += image_data[i+1]
        return merged_image



#
# The I3 grid lazy-load class.
#
# This class will only load the localizations as needed, making
# it quite a bit less memory intensive.
#
class I3GDataLL(I3GData):
    def __init__(self, filename, scale = 4, verbose = True):
        I3GGeneric.__init__(self, 
                            filename,
                            scale = scale,
                            verbose = verbose)

        self.i3_in = readinsight3.I3Reader(filename)
        self.i3data = self.i3_in.nextBlock()
        self.resetFp()

        # Determine film size.
        [image_x, image_y, self.film_l] = getFilmSize(filename, self.i3data)
        self.im_size = [image_x, image_y]

        # Determine what channels the image has.
        self.channels = []
        for i in range(10):
            mask = (self.i3data['c'] == i)
            if mask.sum() > 0:
                self.channels.append(i)

    def close(self):
        self.i3_in.close()

    def dataIsGood(self):
        if(type(self.i3data)==type(numpy.array([]))):
            return True
        else:
            return False

    def getCurrentFrameRange(self):
        return [self.i3data['fr'][0], self.i3data['fr'][-1]]

    def loadDataInFrames(self, fmin = 0, fmax = 500000):
        self.i3data = self.i3_in.getMoleculesInFrameRange(fmin+1, fmax+1)
        self.i3data['fr'] -= 1

    def i3ToXDGridAllChannelsMergedLL(self, grid_fn, verbose):
        self.i3_in.resetFp()
        self.i3data = self.i3_in.nextBlock()
        merged = grid_fn()

        self.i3data = self.i3_in.nextBlock()        
        while(type(self.i3data)==type(numpy.array([]))):
            if verbose:
                sys.stdout.write(".")
                sys.stdout.flush()
            merged += grid_fn()
            self.i3data = self.i3_in.nextBlock()
        if verbose:
            print ""

        return merged

    def i3To2DGridAllChannelsMergedLL(self, zmin = -1000.0, zmax = 1000.0, uncorrected = False, verbose = True):
        def gridFn():
            return self.i3To2DGridAllChannelsMerged(zmin = zmin,
                                                    zmax = zmax,
                                                    uncorrected = uncorrected,
                                                    verbose = verbose)
        return self.i3ToXDGridAllChannelsMergedLL(gridFn, verbose)

    def i3To2DGridSpecificChannelLL(self, channel, zmin = -1000.0, zmax = 1000.0, uncorrected = False, verbose = True):
        def gridFn():
            mask = (self.i3data['c'] == channel)
            self.i3data = i3dtype.maskData(self.i3data, mask)
            return self.i3To2DGridAllChannelsMerged(zmin = zmin,
                                                    zmax = zmax,
                                                    uncorrected = uncorrected,
                                                    verbose = verbose)
        return self.i3ToXDGridAllChannelsMergedLL(gridFn, verbose)

    def i3To3DGridAllChannelsMergedLL(self, z_bins, zmin = -1000.0, zmax = 1000.0, uncorrected = False, verbose = True):
        def gridFn():
            return self.i3To3DGridAllChannelsMerged(z_bins,
                                                    zmin = zmin,
                                                    zmax = zmax,
                                                    uncorrected = uncorrected,
                                                    verbose = verbose)
        return self.i3ToXDGridAllChannelsMergedLL(gridFn, verbose)

    def i3To3DGridSpecificChannelLL(self, channel, z_bins, zmin = -1000.0, zmax = 1000.0, uncorrected = False, verbose = True):
        def gridFn():
            mask = (self.i3data['c'] == channel)
            self.i3data = i3dtype.maskData(self.i3data, mask)
            return self.i3To3DGridAllChannelsMerged(z_bins,
                                                    zmin = zmin,
                                                    zmax = zmax,
                                                    uncorrected = uncorrected,
                                                    verbose = verbose)
        return self.i3ToXDGridAllChannelsMergedLL(gridFn, verbose)

    def nextBlock(self, block_size = 500000):
        self.i3data = self.i3_in.nextBlock(block_size = block_size)
        if(type(self.i3data)==type(numpy.array([]))):
            self.i3data['fr'] -= 1
            return True
        else:
            return False

    def resetFp(self):
        self.i3_in.resetFp()


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
