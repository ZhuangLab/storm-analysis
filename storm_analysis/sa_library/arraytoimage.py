#!/usr/bin/python
#
# Routines for creating images from numpy arrays.
#
# Hazen 7/09
#

from PIL import Image

import colorsys
import numpy

# Save image a 8 bit png with no modifications
def image(data, filename):
    image = numpy.concatenate((data[:,:,None],
                               data[:,:,None],
                               data[:,:,None]),
                              axis = 2)
    image = image.astype(numpy.uint8)
    im = Image.fromarray(image)
    im.save(filename + ".png")

# Creates a image colored by Z from an array of z data.
# Z is assumed normalized to 0.0 - 1.0.
def imageWithZ(zdata, filename, colortable = None, weights = None, channels = None, background = 0):
    if not colortable:
        colortable = zColorTable()
    r = colortable[0].astype(numpy.uint8)
    g = colortable[1].astype(numpy.uint8)
    b = colortable[2].astype(numpy.uint8)

    # create the image
    temp = zdata.copy().astype(numpy.float)
    temp = numpy.round(255.0 * temp)
    img = numpy.zeros((zdata.shape[0], zdata.shape[1], 4))
    img[:,:,0] = r[temp.astype(numpy.uint8)]
    img[:,:,1] = g[temp.astype(numpy.uint8)]
    img[:,:,2] = b[temp.astype(numpy.uint8)]

    # weight the color by the number of counts
    if type(weights) == type(numpy.array([])):
        img[:,:,0] = img[:,:,0] * weights
        img[:,:,1] = img[:,:,1] * weights
        img[:,:,2] = img[:,:,2] * weights

    img = img.astype(numpy.uint8)

    # this sets the image transparency to zero
    img[:,:,3] = 255

    pilimage = Image.fromarray(img, "RGBA")
    pilimage.save(filename + ".png")

# Creates a false-color image of the array data.
# The default color table is gray.
def singleColorImage(data, filename, colortable = None, autoscale = True, xsize = None, ysize = None):

    # default grayscale color table
    r = numpy.arange(0, 256, dtype = numpy.uint8)
    g = numpy.arange(0, 256, dtype = numpy.uint8)
    b = numpy.arange(0, 256, dtype = numpy.uint8)
    if colortable:
        r = colortable[0].astype(numpy.uint8)
        g = colortable[1].astype(numpy.uint8)
        b = colortable[2].astype(numpy.uint8)

    # create the image
    temp = data.copy().astype(numpy.float)
    if autoscale:
        temp = temp - numpy.min(temp)
        if(numpy.max(temp)>0):
            temp = temp/numpy.max(temp)
    temp = numpy.round(255.0 * temp)
    img = numpy.zeros((data.shape[0], data.shape[1], 4), numpy.uint8)
    img[:,:,0] = r[temp.astype(numpy.uint8)]
    img[:,:,1] = g[temp.astype(numpy.uint8)]
    img[:,:,2] = b[temp.astype(numpy.uint8)]

    # this sets the image transparency to zero
    img[:,:,3] = 255

    pilimage = Image.fromarray(img, "RGBA")

    if xsize and ysize:
        pilimage = pilimage.resize([xsize, ysize])
    pilimage.save(filename + ".png")


# Creates a true color image of the array data (a rgb triple).
def trueColorImage(data, filename):
    # create the image
    sx = 0
    sy = 0
    rgb = []
    for datum in data:
        if type(datum) == type(numpy.array([])):
            sx = datum.shape[0]
            sy = datum.shape[1]
            temp = datum.copy().astype(numpy.float)
            temp = numpy.round(255.0 * temp)
            rgb.append(temp)
        else:
            rgb.append(0)

    img = numpy.zeros((sx, sy, 4), numpy.uint8)
    if type(rgb[0]) == type(numpy.array([])):
        img[:,:,0] = rgb[0].astype(numpy.uint8)
    if type(rgb[1]) == type(numpy.array([])):
        img[:,:,1] = rgb[1].astype(numpy.uint8)
    if type(rgb[2]) == type(numpy.array([])):
        img[:,:,2] = rgb[2].astype(numpy.uint8)

    # white background (sort of)
    # mask = (numpy.sum(img, axis = 2) == 0)
    # img[mask,:] = 255

    # this sets the image transparency to zero
    img[:,:,3] = 255

    pilimage = Image.fromarray(img, "RGBA")
    pilimage.save(filename + ".png")


# Default Z color table
def zColorTable():
    h = 1.0 - (1.0/255.0) * numpy.arange(0.0, 256.0)
    l = 0.5 * numpy.ones(256)
    s = 1.0 * numpy.ones(256)
    r = numpy.ones(256, dtype = numpy.float32)
    g = numpy.ones(256, dtype = numpy.float32)
    b = numpy.ones(256, dtype = numpy.float32)
    for i in range(256):
        [_r, _g, _b] = colorsys.hls_to_rgb(h[i], l[i], s[i])
        r[i] = 255.0 * _r
        g[i] = 255.0 * _g
        b[i] = 255.0 * _b
    return [r, g, b]


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
