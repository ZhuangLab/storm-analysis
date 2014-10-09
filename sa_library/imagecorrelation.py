#!/usr/local/bin/python
#
# Classes and functions for image correlation.
#
# Hazen 07/09
#

import numpy
import scipy
import scipy.signal

import gaussfit

def absIntRound(num):
    return abs(int(round(num)))

def crop2DImages(image1, image2, dx, dy):
    [size_x, size_y] = image1.shape
    adx = absIntRound(dx)
    ady = absIntRound(dy)
    if dx >= 0 and dy >= 0:
        image1 = image1[:size_x-adx,:size_y-ady]
        image2 = image2[adx:       ,ady:       ]
    if dx >= 0 and dy < 0:
        image1 = image1[:size_x-adx,ady:       ]
        image2 = image2[adx:       ,:size_y-ady]
    if dx < 0 and dy >= 0:
        image1 = image1[adx:       ,:size_y-ady]
        image2 = image2[:size_x-adx,ady:       ]
    if dx < 0 and dy < 0:
        image1 = image1[adx:       ,ady:       ]
        image2 = image2[:size_x-adx,:size_y-ady]
    return [image1, image2]

def crop3DImages(image1, image2, dx, dy):
    [size_x, size_y, size_z] = image1.shape
    adx = absIntRound(dx)
    ady = absIntRound(dy)
    if dx >= 0 and dy >= 0:
        image1 = image1[:size_x-adx,:size_y-ady,:]
        image2 = image2[adx:       ,ady:       ,:]
    if dx >= 0 and dy < 0:
        image1 = image1[:size_x-adx,ady:       ,:]
        image2 = image2[adx:       ,:size_y-ady,:]
    if dx < 0 and dy >= 0:
        image1 = image1[adx:       ,:size_y-ady,:]
        image2 = image2[:size_x-adx,ady:       ,:]
    if dx < 0 and dy < 0:
        image1 = image1[adx:       ,ady:       ,:]
        image2 = image2[:size_x-adx,:size_y-ady,:]
    return [image1, image2]

def xyCorrelate(image1, image2):
    return scipy.signal.fftconvolve(image1, image2[::-1, ::-1], mode="same")

# Note that the search is limited to a X by X region
# in the center of the overlap between the two images.
def xyOffset(image1, image2, scale, center = None):
    image1 = image1 - numpy.median(image1)
    image2 = image2 - numpy.median(image2)

    if 0:
        import arraytoimage
        print "1:", numpy.max(image1), "2:", numpy.max(image2)
        arraytoimage.singleColorImage(image1, "corr_image1")
        arraytoimage.singleColorImage(image2, "corr_image2")

    result = xyCorrelate(image1, image2)
    mx = int(round(0.5 * result.shape[0]))
    my = int(round(0.5 * result.shape[1]))
    s_size = 30 * int(scale)
    sx = s_size
    sy = s_size
    if mx < (s_size + 5):
        sx = mx - 5
    if my < (s_size + 5):
        sy = my - 5
    if type(center) == type([]):
        rx = round(center[0]) + sx
        ry = round(center[1]) + sy
        if 0:
            print rx
            print numpy.argmax(numpy.max(result[mx-sx:mx+sx+1,my-sy:my+sy+1], axis = 1))
            print ry
            print numpy.argmax(numpy.max(result[mx-sx:mx+sx+1,my-sy:my+sy+1], axis = 0))
    else:
        rx = numpy.argmax(numpy.max(result[mx-sx:mx+sx+1,my-sy:my+sy+1], axis = 1))
        ry = numpy.argmax(numpy.max(result[mx-sx:mx+sx+1,my-sy:my+sy+1], axis = 0))
    rx += (mx - sx)
    ry += (my - sy)
    [fit, success] = gaussfit.fitSymmetricGaussian(result[rx-5:rx+6,ry-5:ry+6], 2.0)

    if 0:
        fit_fn = gaussfit.symmetricGaussian(*fit)
        surf = numpy.zeros((11,11))
        fp = open("fit.txt", "w")
        fp.write("x, y, correlation, fit\n")
        for x in range(11):
            for y in range(11):
                surf[x,y] = fit_fn(x,y)
                fp.write(str(x) + ", " + str(y) + ", " + str(result[rx-5+x,ry-5+y]) + ", " + str(surf[x,y]) + "\n")
        fp.close()

    if (fit[0] == fit[1]) or (fit[2] < 2.0) or (fit[2] > 8.0) or (fit[3] < 2.0) or (fit[3] > 8.0):
        print "Bad fit center:", fit[2], fit[3]
        success = False
    return [result[mx-sx:mx+sx+1,my-sy:my+sy+1],
            rx + fit[2] - 5.0 - 0.5 * result.shape[0],
            ry + fit[3] - 5.0 - 0.5 * result.shape[1],
            success]

def xyOffsetWithDxDy(image1, image2, dx, dy, scale):
    [cimage1, cimage2] = crop2DImages(image1, image2, dx, dy)
    return xyOffset(cimage1, cimage2, scale)

def zOffset(image1, image2):
    image1 = image1 - numpy.median(image1)
    image2 = image2 - numpy.median(image2)
    [size_x, size_y, size_z] = image1.shape
    corr = numpy.zeros(2*size_z-1)
    for i in range(size_z-1):
        corr[i] = numpy.sum(image1[:,:,size_z-1-i:size_z] * image2[:,:,:i+1])/float(i+1)
    for i in range(size_z):
        corr[size_z-1+i] = numpy.sum(image1[:,:,:size_z-i] * image2[:,:,i:size_z])/float(size_z-i)
    
    # this handles data that is actually 2D
    number_non_zero = 0
    which_non_zero = 0
    for i in range(2*size_z-1):
        if (corr[i] > 0.0):
            which_non_zero = i
            number_non_zero += 1

    if (number_non_zero > 1):
        [fit, success] = gaussfit.fitLorentzian(corr)
        lorentzian_fit = gaussfit.lorentzian(*fit)(*numpy.indices(corr.shape))
    else:
        fit = [0, 0, which_non_zero]
        success = True
        lorentzian_fit = corr

    return [corr, lorentzian_fit, fit[2] - size_z + 1, success]

def zOffsetWithDxDy(image1, image2, dx, dy):
    [cimage1, cimage2] = crop3DImages(image1, image2, dx, dy)
    return zOffset(cimage1, cimage2)



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
