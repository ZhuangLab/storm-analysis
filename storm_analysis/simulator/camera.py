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

import storm_analysis.simulator.simbase as simbase

class Camera(simbase.SimBase):
    """
    Converts the image from photons to counts and adds camera
    noise, baseline, etc.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, baseline):
        simbase.SimBase.__init__(self, sim_fp, x_size, y_size, i3_data)
        self.baseline = baseline


class Ideal(Camera):
    """
    Perfect camera with only shot noise.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, baseline):
        Camera.__init__(self, sim_fp, x_size, y_size, i3_data, baseline)
        self.saveJSON({"camera" : {"class" : "Ideal",
                                   "baseline" : str(baseline)}})

    def readImage(self, image):
        return numpy.random.poisson(image) + self.baseline


class EMCCD(Camera):
    """
    A camera with EMCCD gain.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, baseline, emccd_gain = 30.0, preamp_gain = 1.0/5.0, read_noise = 20.0):
        Camera.__init__(self, sim_fp, x_size, y_size, i3_data, baseline)
        self.emccd_gain = emccd_gain
        self.preamp_gain = preamp_gain
        self.read_noise = read_noise
        self.saveJSON({"camera" : {"class" : "Ideal",
                                   "baseline" : str(baseline),
                                   "emccd_gain" : str(emccd_gain),
                                   "preamp_gain" : str(preamp_gain),
                                   "read_noise" : str(read_noise)}})

    def readImage(self, image):
        """
        I think this is right..
        """

        # Detected image.
        image = numpy.random.poisson(image).astype(numpy.float64)
        
        # EMCCD part.
        image = numpy.random.poisson(self.emccd_gain * image).astype(numpy.float64)

        # Add read-noise.
        image += numpy.random.normal(scale = self.read_noise, size = image.shape)

        # Pre-amp.
        image = image * self.preamp_gain

        return image


class SCMOS(Camera):
    """
    A sCMOS camera. The sCMOS calibration data needs to be the same size as
    the simulated images.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, baseline, scmos_cal):
        Camera.__init__(self, sim_fp, x_size, y_size, i3_data, baseline)
        [self.offset, variance, self.gain] = numpy.load(scmos_cal)
        self.std_dev = numpy.sqrt(variance)

        if (self.offset.shape[0] != x_size) or (self.offset.shape[1] != y_size):
            raise simbase.SimException("sCMOS calibration data size does not match the image size.")
        
        self.saveJSON({"camera" : {"class" : "Ideal",
                                   "baseline" : str(baseline),
                                   "scmos_cal" : str(emccd_gain)}})

    def readImage(self, image):

        # Detected image.
        image = numpy.random.poisson(image)

        # Multiply by pixel dependent gain.
        image = self.gain * image

        # Add pixel dependent noise.
        image += numpy.random.normal(scale = self.std_dev)

        # Add pixel dependent offset.'
        image += self.offset

        return image

        
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
