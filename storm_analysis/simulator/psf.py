#!/usr/bin/python
#
# Classes for creating different kinds of PSFs.
#
# Hazen 11/16
#

import numpy
import random

import storm_analysis.simulator.draw_gaussians_c as dg
import storm_analysis.simulator.pupil_math as pupilMath
import storm_analysis.simulator.simbase as simbase


class PSF(simbase.SimBase):
    """
    Draws the emitter PSFs on an image.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, nm_per_pixel):
        simbase.SimBase.__init__(self, sim_fp, x_size, y_size, i3_data)
        self.nm_per_pixel = nm_per_pixel


class GaussianPSF(PSF):
    """
    Gaussian PSF.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, nm_per_pixel):
        PSF.__init__(self, sim_fp, x_size, y_size, i3_data, nm_per_pixel)
        self.saveJSON({"psf" : {"class" : "GaussianPSF",
                                "nm_per_pixel" : str(nm_per_pixel)}})

    def getPSFs(self, i3_data):
        x = i3_data['x'] - 1.0
        y = i3_data['y'] - 1.0
        h = i3_data['h']

        ax = i3_data['ax']
        w = i3_data['w']
        sx = 0.5*numpy.sqrt(w*w/ax)/self.nm_per_pixel
        sy = 0.5*numpy.sqrt(w*w*ax)/self.nm_per_pixel

        return dg.drawGaussians((self.x_size, self.y_size),
                                numpy.concatenate((x[:,None],
                                                   y[:,None],
                                                   h[:,None],
                                                   sx[:,None],
                                                   sy[:,None]),
                                                  axis = 1))

class PupilFunction(PSF):
    """
    PSF using the pupil function approach.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, nm_per_pixel, zmn, wavelength = 600, refractive_index = 1.5, numerical_aperture = 1.4):
        """
        zmn is a list of lists containing the zernike mode terms, e.g.
            [[1.3, 2, 2]] for pure astigmatism.
        wavelength is the mean emission wavelength in nm.
        """
        PSF.__init__(self, sim_fp, x_size, y_size, i3_data, nm_per_pixel)
        self.saveJSON({"psf" : {"class" : "PupilFunction",
                                "nm_per_pixel" : str(nm_per_pixel),
                                "numerical_aperture" : str(numerical_aperture),
                                "refactrive_index" : str(refractive_index),
                                "wavelength" : str(wavelength),
                                "zmn" : str(zmn)}})

        self.geo = pupilMath.Geometry(int(4.0/(nm_per_pixel * 0.001)),
                                      nm_per_pixel * 0.001,
                                      wavelength * 0.001,
                                      refractive_index,
                                      numerical_aperture)
        self.pf = self.geo.createFromZernike(1.0, zmn)
        self.psf_size = self.geo.r.shape[0]
        self.margin = int(self.psf_size/2) + 1

        self.im_size_x = self.x_size + 2 * self.margin
        self.im_size_y = self.y_size + 2 * self.margin

    def getPSFs(self, i3_data):
        image = numpy.zeros((self.im_size_x, self.im_size_y))
        x = i3_data['x']         # Pixels
        y = i3_data['y']         # Pixels
        z = i3_data['z']*0.001   # Expected to be in nanometers.
        h = i3_data['h']

        dx = x - numpy.floor(x)
        dy = y - numpy.floor(y)

        for i in range(x.size):

            ix = int(x[i])
            iy = int(y[i])
            
            if (ix >= 0.0) and (ix < self.x_size) and (iy >= 0.0) and (iy < self.y_size):

                # Shift to the desired z value.
                defocused = self.geo.changeFocus(self.pf, z[i])

                # Translate to the correct sub-pixel position.
                #translated = defocused
                translated = self.geo.translatePf(defocused, dx[i], dy[i])

                # Get real-space intensity.
                psf = pupilMath.intensity(pupilMath.toRealSpace(translated))
                
                image[ix:ix+self.psf_size,iy:iy+self.psf_size] += h[i] * psf

        return image[self.margin:self.margin+self.x_size,self.margin:self.margin+self.y_size]


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
