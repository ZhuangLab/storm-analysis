#!/usr/bin/env python
"""
Classes simulating different kinds of cameras.

These are also responsible for adding the Poisson noise
to the image (if desired).

Hazen 11/16
"""

import numpy
import pickle
import random

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.simulator.simbase as simbase

class Camera(simbase.SimBase):
    """
    Converts the image from photons to counts and adds camera
    noise, baseline, etc.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data):
        simbase.SimBase.__init__(self, sim_fp, x_size, y_size, i3_data)


class Ideal(Camera):
    """
    Perfect camera with only shot noise.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, baseline, gain = 1.0):
        super(Ideal, self).__init__(sim_fp, x_size, y_size, i3_data)
        
        self.baseline = baseline
        self.gain = gain
        self.saveJSON({"camera" : {"class" : "Ideal",
                                   "baseline" : str(baseline),
                                   "gain" : str(gain)}})

    def readImage(self, image):
        return self.gain*numpy.random.poisson(image) + self.baseline


class IdealNoNoise(Camera):
    """
    Perfect camera without shot noise.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, baseline, gain = 1.0):
        super(IdealNoNoise, self).__init__(sim_fp, x_size, y_size, i3_data)
        
        self.baseline = baseline
        self.gain = gain
        self.saveJSON({"camera" : {"class" : "IdealNoNoise",
                                   "baseline" : str(baseline),
                                   "gain" : str(gain)}})

    def readImage(self, image):
        return self.gain*image + self.baseline

    
class EMCCD(Camera):
    """
    A camera with EMCCD gain.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, baseline, emccd_gain = 30.0, preamp_gain = 1.0/5.0, read_noise = 20.0):
        super(EMCCD, self).__init__(sim_fp, x_size, y_size, i3_data)
        
        self.baseline = baseline
        self.emccd_gain = emccd_gain
        self.preamp_gain = preamp_gain
        self.read_noise = read_noise
        self.saveJSON({"camera" : {"class" : "EMCCD",
                                   "baseline" : str(baseline),
                                   "emccd_gain" : str(emccd_gain),
                                   "preamp_gain" : str(preamp_gain),
                                   "read_noise" : str(read_noise)}})

    def readImage(self, image):

        # Detected image.
        image = numpy.random.poisson(image)
        
        # EMCCD part.
        for i in range(image.shape[0]):
            for j in range(image.shape[1]):
                image[i,j] = numpy.sum(numpy.random.exponential(self.emccd_gain, image[i,j]))
        image = image.astype(numpy.float64)
        
        # Add read-noise.
        image += numpy.random.normal(scale = self.read_noise, size = image.shape)

        # Pre-amp.
        image = image * self.preamp_gain

        return image + self.baseline


class SCMOS(Camera):
    """
    A sCMOS camera. The sCMOS calibration data needs to be the same size as
    the simulated images.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, scmos_cal):
        super(SCMOS, self).__init__(sim_fp, x_size, y_size, i3_data)

        # We transpose the calibration data here because the images we process
        # are the transpose of the images that are saved for the analysis.
        #
        [self.offset, variance, self.gain, self.rqe] = map(numpy.transpose,
                                                           analysisIO.loadCMOSCalibration(scmos_cal))
        self.std_dev = numpy.sqrt(variance)

        if (self.offset.shape[0] != x_size) or (self.offset.shape[1] != y_size):
            raise simbase.SimException("sCMOS calibration data size does not match the image size.")
        
        self.saveJSON({"camera" : {"class" : "SCMOS",
                                   "scmos_cal" : scmos_cal}})

    def readImage(self, image):

        # Multiply by relative QE.
        image = image * self.rqe
        
        # Detected image.
        image = numpy.random.poisson(image)

        # Multiply by pixel dependent gain.
        image = self.gain * image

        # Add pixel dependent noise.
        image += numpy.random.normal(scale = self.std_dev)

        # Add pixel dependent offset.
        image += self.offset
        
        return image


def createSCMOSCalibration(x_size, y_size, gain, read_noise, hot_fraction = 0.05, hot_lambda = None, offset = 100.0):
    """
    Create a simulated calibration file for a sCMOS camera.

    The camera is modeled as having a gaussian distribution of read noise values
    with mean = 'read_noise' and sigma = 0.1 * read_noise. In addition a small
    fraction of the pixels are 'hot' and have an additional exponential distributed
    value added to their read noise.

    x_size - Size in x in pixels.
    y_size - Size in y in pixels.
    gain - Camera gain in ADU/electrons.
    read_noise - Read noise in electrons.
    hot_fraction - The fraction of the pixels that are 'hot'.
    hot_lambda - The lambda for the exponential distribution of the hot pixels.
    offset - Camera baseline in ADU.

    returns [offset, variance, gain, rqe]
    """
    offset = numpy.zeros((x_size, y_size)) + offset
    gain = gain * numpy.ones((x_size, y_size))
    rqe = numpy.ones((x_size, y_size))

    # We're squaring everything as calibration file contains the variance, not the
    # standard deviation.
    #
    r_noise = numpy.random.normal(loc = read_noise, scale = 0.1 * read_noise, size = (x_size, y_size))
    r_noise = r_noise * r_noise

    # Use 10x the read noise if hot_lambda is not specified.
    if hot_lambda is None:
        hot_lambda = 10.0 * read_noise
    
    e_noise = numpy.random.exponential(scale = hot_lambda, size = (x_size, y_size))
    e_noise = e_noise * e_noise

    mask = (numpy.random.uniform(size = (x_size, y_size)) > hot_fraction)
    e_noise[mask] = 0.0

    variance = r_noise + e_noise

    # Convert variance to ADU^2
    variance = variance * gain * gain

    return [offset, variance, gain, rqe]
    

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
