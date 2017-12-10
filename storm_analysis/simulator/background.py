#!/usr/bin/env python
"""
Classes for creating different kinds of backgrounds.

Hazen 11/16
"""

import numpy
import random

import storm_analysis.simulator.draw_gaussians_c as dg
import storm_analysis.simulator.simbase as simbase


class Background(simbase.SimBase):
    """
    Generate a background image (in photons).
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data):
        simbase.SimBase.__init__(self, sim_fp, x_size, y_size, i3_data)

    def getBackground(self, frame):
        return self.bg_image

    def getEmitterBackground(self, i3_data_in):
        i3_data = numpy.copy(i3_data_in)
        for i in range(i3_data['x'].size):
            x = int(round(i3_data['x'][i]))
            y = int(round(i3_data['y'][i]))
            i3_data['bg'][i] = self.bg_image[x,y]
        return i3_data

    
class GaussianBackground(Background):

    def __init__(self, sim_fp, x_size, y_size, i3_data, photons = 100, sigma = 100.0):
        Background.__init__(self, sim_fp, x_size, y_size, i3_data)
        self.saveJSON({"background" : {"class" : "GaussianBackground",
                                       "photons" : str(photons),
                                       "sigma" : str(sigma)}})
        self.bg_image = dg.drawGaussiansXY((x_size, y_size),
                                           numpy.array([0.5*x_size]),
                                           numpy.array([0.5*y_size]),
                                           sigma = sigma)
        self.bg_image = photons * self.bg_image/numpy.max(self.bg_image)


class SlopedBackground(Background):

    def __init__(self, sim_fp, x_size, y_size, i3_data, slope = 0.1, offset = 0.0):
        Background.__init__(self, sim_fp, x_size, y_size, i3_data)
        self.saveJSON({"background" : {"class" : "SlopedBackground",
                                       "offset" : str(offset),
                                       "slope" : str(slope)}})
        self.bg_image = numpy.zeros((x_size, y_size)) + offset

        slope_arr = numpy.arange(x_size) * slope
        for i in range(y_size):
            self.bg_image[:,i] += slope_arr
            
    
class UniformBackground(Background):

    def __init__(self, sim_fp, x_size, y_size, i3_data, photons = 100):
        Background.__init__(self, sim_fp, x_size, y_size, i3_data)
        self.saveJSON({"background" : {"class" : "UniformBackground",
                                       "photons" : str(photons)}})
        self.bg_image = numpy.ones((x_size, y_size)) * photons



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
