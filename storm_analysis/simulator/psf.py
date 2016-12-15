#!/usr/bin/python
#
# Classes for creating different kinds of PSFs.
#
# Hazen 11/16
#

import numpy
import pickle
import random

import storm_analysis.spliner.spline_to_psf as splineToPSF

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

        if ((self.psf_size%2)==0):
            self.margin = int(self.psf_size/2) + 1
        else:
            self.margin = int(self.psf_size/2) + 2

        print("psf size", self.psf_size)

        self.im_size_x = self.x_size + 2 * self.margin
        self.im_size_y = self.y_size + 2 * self.margin

    def getPSFs(self, i3_data):
        """
        The expected form for the i3 data fields are x,y in pixels and z in nanometers.
        """
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


class Spline2D(splineToPSF.SplineToPSF2D):
    """
    2D spline with non-zero offsets in x, y.
    """
    def __init__(self, spline_file):
        splineToPSF.SplineToPSF2D.__init__(self, spline_file)
        self.psf_size = int((self.spline_size - 1)/2) - 1
        
    def getPSF(self, z_value, dx, dy):
        """
        z_value is ignored, it is only a parameter so that
        the signature matches Spline3D.

        dx, dy are in the range 0.0 - 1.0.
        """
        psf = numpy.zeros((self.psf_size, self.psf_size))
        if(((self.psf_size+1)%2) == 0):
            for x in range(self.psf_size):
                for y in range(self.psf_size):
                    psf[y,x] = self.spline.f(float(2*(y+dy)),
                                             float(2*(x+dx)))
        else:
            for x in range(self.psf_size):
                for y in range(self.psf_size):
                    psf[y,x] = self.spline.f(float(2*(y+dy)) + 1.0,
                                             float(2*(x+dx)) + 1.0)
            
        return psf


class Spline3D(splineToPSF.SplineToPSF3D):
    """
    3D spline with non-zero offsets in x, y.
    """
    def __init__(self, spline_file):
        splineToPSF.SplineToPSF3D.__init__(self, spline_file)
        self.psf_size = int((self.spline_size - 1)/2) - 1
        
    def getPSF(self, z_value, dx, dy):
        """
        z_value needs to be inside the z range covered by the spline.

        dx, dy are in the range 0.0 - 1.0.
        """
        scaled_z = self.getScaledZ(z_value)
                
        psf = numpy.zeros((self.psf_size, self.psf_size))
        if(((self.psf_size+1)%2) == 0):
            for x in range(self.psf_size):
                for y in range(self.psf_size):
                    psf[y,x] = self.spline.f(scaled_z,
                                             float(2*(y+dy)),
                                             float(2*(x+dx)))
        else:
            for x in range(self.psf_size):
                for y in range(self.psf_size):
                    psf[y,x] = self.spline.f(scaled_z,
                                             float(2*(y+dy)) + 1.0,
                                             float(2*(x+dx)) + 1.0)
            
        return psf

    
class Spline(PSF):
    """
    PSF from a (cubic) spline.
    """
    def __init__(self, sim_fp, x_size, y_size, i3_data, nm_per_pixel, spline_file):
        """
        spline_file is the name of a .spline file as generated by spliner/psf_to_spline.py.

        Splines are always pixel based, so the spline pixel size should match what you
        want for the simulation. The pixel size of the spline is recorded in the spline
        file and in theory that could be used to adjust appropriately depending on the
        desired pixel size for the simulation, but that is not currently being done.
        """
        PSF.__init__(self, sim_fp, x_size, y_size, i3_data, nm_per_pixel)
        self.saveJSON({"psf" : {"class" : "Spline",
                                "spline_file" : spline_file}})
        spline_data = pickle.load(open(spline_file, 'rb'))
        if (spline_data["type"] == "3D"):
            self.spline = Spline3D(spline_data)
        else:
            self.spline = Spline2D(spline_data)

        self.psf_size = self.spline.psf_size

        if ((self.psf_size%2)==0):
            self.margin = int(self.psf_size/2)
        else:
            self.margin = int(self.psf_size/2) + 1

        self.im_size_x = self.x_size + 2*self.margin
        self.im_size_y = self.y_size + 2*self.margin

    def getPSFs(self, i3_data):
        """
        The expected form for the i3 data fields are x,y in pixels and z in nanometers.
        """
        image = numpy.zeros((self.im_size_x, self.im_size_y))
        x = i3_data['x']   # Pixels
        y = i3_data['y']   # Pixels
        z = i3_data['z']   # Expected to be in nanometers.
        h = i3_data['h']

        dx = 1.0 - (x - numpy.floor(x))
        dy = 1.0 - (y - numpy.floor(y))

        for i in range(x.size):

            ix = int(x[i])
            iy = int(y[i])
            
            if (ix >= 0.0) and (ix < self.x_size) and (iy >= 0.0) and (iy < self.y_size):

                # Calculate psf, dx and dy are transposed to match the spline coordinate system.
                psf = self.spline.getPSF(z[i], dy[i], dx[i])

                # Scale to correct height.
                psf = h[i] * psf/numpy.sum(psf)

                # Add to image.
                image[ix:ix+self.psf_size,iy:iy+self.psf_size] += psf

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
