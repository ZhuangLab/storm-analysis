#!/usr/bin/env python
#
# 02/16
#
# Simple Python interface to rolling_ball_lib.c.
#
# Hazen
#

import numpy
import os
import scipy
import scipy.ndimage

import ctypes
import numpy
from numpy.ctypeslib import ndpointer

import sa_library.loadclib as loadclib

rball = loadclib.loadCLibrary(os.path.dirname(__file__), "rolling_ball_lib")

# C interface definition
rball.estimateBg.argtypes = [ndpointer(dtype=numpy.float64),
                             ndpointer(dtype=numpy.float64),
                             ctypes.c_int,
                             ctypes.c_int]
rball.init.argtypes = [ndpointer(dtype=numpy.float64), 
                       ctypes.c_int]

class CRollingBall(object):

    def __init__(self, ball_radius, smoothing_sigma):
        self.ball_radius = ball_radius
        self.smoothing_sigma = smoothing_sigma

        ball_size = int(round(ball_radius * 0.5))
        br = ball_radius * ball_radius
        ball = numpy.zeros((2*ball_size+1, 2*ball_size+1))
        for x in range(2*ball_size+1):
            dx = x - ball_size
            for y in range(2*ball_size+1):
                dy = y - ball_size
                ball[x,y] = br - (dx * dx + dy * dy)
        ball = numpy.sqrt(ball)

        rball.init(ball, 2*ball_size+1)

    def estimateBG(self, image):
        image = image.astype(numpy.float)
        sm_image = scipy.ndimage.filters.gaussian_filter(image, self.smoothing_sigma)
        ball_image = numpy.zeros((image.shape))
        rball.estimateBg(sm_image, ball_image, sm_image.shape[0], sm_image.shape[1])
        return ball_image + self.ball_radius

    def removeBG(self, image):
        return image - self.estimateBG(image)


if (__name__ == "__main__"):

    # A very simple test.
    test = 100.0 * numpy.ones((100,100))
    rb = CRollingBall(10, 0.5)
    test = rb.removeBG(test)
    print numpy.min(test), numpy.max(test)
    

#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
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
