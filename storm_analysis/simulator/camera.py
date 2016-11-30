#!/usr/bin/python
#
# Classes simulating different kinds of cameras.
#
# These are also responsible for adding the Poisson noise
# to the image (if desired).
#
# Hazen 11/16
#

import numpy
import random

import storm_analysis.simulator.draw_gaussians_c as dg


class Camera(object):
    """
    Converts the image from photons to counts and adds camera
    noise, baseline, etc.
    """
    def __init__(self, x_size, y_size, baseline):
        self.baseline = baseline
        self.x_size = x_size
        self.y_size = y_size


class Ideal(Camera):
    """
    Perfect camera with only shot noise.
    """
    def __init__(self, x_size, y_size, baseline):
        Camera.__init__(self, x_size, y_size)

    def readImage(self, image):
        return numpy.random.poisson(image) + self.baseline

        
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
