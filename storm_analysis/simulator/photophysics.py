#!/usr/bin/python
#
# Classes for simulating different dye photophysics.
#
# Hazen 11/16
#

import numpy
import random


class PhotoPhysics(object):
    """
    Returns location and intensity (peak height in photons) of
    the emitters that are on in the current frame.
    """
    def __init__(self, x_size, y_size, i3_data):
        self.i3_data = i3_data
        self.x_size = x_size
        self.y_size = y_size

class AlwaysOn(PhotoPhysics):
    """
    All the emitters are on all the time.
    """
    def __init__(self, x_size, y_size, i3_data, intensity = 100):
        PhotoPhysics.__init__(self, x_size, y_size, i3_data)
        self.i3_data['h'] = intensity

    def getEmitters(self, frame):
        return self.i3_data


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
