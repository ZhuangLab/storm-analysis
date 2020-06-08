#!/usr/bin/env python
"""
Classes for creating different kinds of PSFs.

Note that as a side-effect the getPSFs() method is expected
to also set the correct value of the 'h' field based on the
PSF and the 'a' field.

Hazen 11/16
"""

import numpy
import pickle
import random
import scipy
import scipy.fftpack

import storm_analysis.spliner.spline_to_psf as splineToPSF

import storm_analysis.simulator.draw_gaussians_c as dg
import storm_analysis.simulator.pupil_math as pupilMath
import storm_analysis.simulator.simbase as simbase

# These are here to make it easier to create pupil functions that match
# those used in the simulations.
pf_size = 128
pf_wavelength = 600
pf_refractive_index = 1.5
pf_numerical_aperture = 1.4

class PSF(simbase.SimBase):
    """
    Draws the emitter PSFs on an image.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, nm_per_pixel):
        super(PSF, self).__init__(sim_fp, x_size, y_size, h5_data)
        self.nm_per_pixel = nm_per_pixel


class DHPSF(PSF):
    """
    A very simplistic approximation of the double helix PSF.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, nm_per_pixel, z_range = 0.75):
        """
        z_range - (half) Z range in microns.
        """
        super(DHPSF, self).__init__(sim_fp, x_size, y_size, h5_data, nm_per_pixel)
        self.z_max = z_range
        self.z_min = -z_range

        self.saveJSON({"psf" : {"class" : "DHPSF",
                                "z_range" : str(z_range)}})

    def getPSFs(self, h5_data):
        x = h5_data['x']
        y = h5_data['y']
        a = h5_data['sum']

        sx = h5_data['ysigma']
        sy = h5_data['xsigma']

        h = 0.5*a/(2.0 * numpy.pi * sx * sy)
        h5_data['height'] = h

        angle = numpy.pi * 0.9 * ((h5_data['z'] - self.z_min)/(self.z_max - self.z_min))
        dx = 2.0 * sx * numpy.cos(angle)
        dy = 2.0 * sy * numpy.sin(angle)

        x1 = x + dx
        y1 = y + dy
        x2 = x - dx
        y2 = y - dy

        x = numpy.concatenate((x1, x2))
        y = numpy.concatenate((y1, y2))
        h = numpy.concatenate((h, h))
        sx = numpy.concatenate((sx, sx))
        sy = numpy.concatenate((sy, sy))
        
        return dg.drawGaussians((self.x_size, self.y_size),
                                numpy.stack((x, y, h, sx, sy), axis = 1),
                                res = 5)

    
class GaussianPSF(PSF):
    """
    Gaussian PSF.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, nm_per_pixel):
        super(GaussianPSF, self).__init__(sim_fp, x_size, y_size, h5_data, nm_per_pixel)
        self.saveJSON({"psf" : {"class" : "GaussianPSF",
                                "nm_per_pixel" : str(nm_per_pixel)}})

    def getPSFs(self, h5_data):
        x = h5_data['x']
        y = h5_data['y']
        a = h5_data['sum']

        sx = h5_data['ysigma']
        sy = h5_data['xsigma']

        h = a/(2.0 * numpy.pi * sx * sy)
        h5_data['height'] = h

        return dg.drawGaussians((self.x_size, self.y_size),
                                numpy.concatenate((x[:,None],
                                                   y[:,None],
                                                   h[:,None],
                                                   sx[:,None],
                                                   sy[:,None]),
                                                  axis = 1),
                                res = 5)


class PupilFunctionBase(PSF):
    """
    PSF using the pupil function approach.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, nm_per_pixel, pf_size):
        super(PupilFunctionBase, self).__init__(sim_fp, x_size, y_size, h5_data, nm_per_pixel)

        self.otf_scaler = None
        self.psf_size = pf_size
        self.sample_index = None
        self.z_center = None

        self.margin = int(self.psf_size/2) + 1

        self.im_size_x = self.x_size + 2 * self.margin
        self.im_size_y = self.y_size + 2 * self.margin

    def getAberrationFn(self, z_value):
        return None

    def getPSFs(self, h5_data):
        """
        The expected form for the h5 data fields are x,y in pixels and z in microns.
        """
        image = numpy.zeros((self.im_size_x, self.im_size_y))
        x = h5_data['x']+1       # Pixels
        y = h5_data['y']+1       # Pixels
        z = h5_data['z']         # Expected to be in microns.
        a = h5_data['sum']

        h5_data['height'] = numpy.zeros(a.size)

        dx = x - numpy.floor(x)
        dy = y - numpy.floor(y)

        for i in range(x.size):

            ix = int(x[i])
            iy = int(y[i])
            
            if (ix >= 0.0) and (ix < self.x_size) and (iy >= 0.0) and (iy < self.y_size):

                pf = self.pf
                
                # Apply aberration function if available.
                ab_fn = self.getAberrationFn(z[i])                
                if ab_fn is not None:
                    pf = pf*ab_fn

                # Shift to the desired z value.
                defocused = self.geo.changeFocus(pf, z[i])

                # Translate to the correct sub-pixel position.
                #translated = defocused
                translated = self.geo.translatePf(defocused, dx[i], dy[i])

                # Get real-space intensity.
                psf = pupilMath.intensity(pupilMath.toRealSpace(translated)) * a[i]

                # Apply OTF scaling if requested.
                if self.otf_scaler is not None:
                    otf = scipy.fftpack.fftshift(scipy.fftpack.fft2(psf))
                    otf_scaled = otf * self.otf_scaler
                    psf = numpy.abs(scipy.fftpack.ifft2(otf_scaled))
                
                h5_data['height'][i] = numpy.max(psf)

                image[ix:ix+self.psf_size,iy:iy+self.psf_size] += psf

        return image[self.margin:self.margin+self.x_size,self.margin:self.margin+self.y_size]
        
    
class PupilFunction(PupilFunctionBase):
    """
    PSF using the pupil function approach. This is the most basic version.

    Note that by default this uses pupilMath.GeometrySim class which (more 
    or less) emulates OTF scaling by using a pixel size that is half the 
    specified pixel size to create the PSF.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, nm_per_pixel, zmn,
                 wavelength = pf_wavelength,
                 refractive_index = pf_refractive_index,
                 numerical_aperture = pf_numerical_aperture,
                 pf_size = pf_size,
                 geo_sim_pf = True):
        """
        zmn is a list of lists containing the zernike mode terms, e.g.
            [[1.3, 2, 2]] for pure astigmatism.
        wavelength is the mean emission wavelength in nm.
        """
        super(PupilFunction, self).__init__(sim_fp, x_size, y_size, h5_data, nm_per_pixel, pf_size)
        self.saveJSON({"psf" : {"class" : "PupilFunction",
                                "geo_sim_pf" : str(geo_sim_pf),
                                "nm_per_pixel" : str(nm_per_pixel),
                                "numerical_aperture" : str(numerical_aperture),
                                "pf_size" : str(pf_size),
                                "refactrive_index" : str(refractive_index),
                                "wavelength" : str(wavelength),
                                "zmn" : str(zmn)}})

        if geo_sim_pf:
            self.geo = pupilMath.GeometrySim(pf_size,
                                             nm_per_pixel * 0.001,
                                             wavelength * 0.001,
                                             refractive_index,
                                             numerical_aperture)
        else:
            self.geo = pupilMath.Geometry(pf_size,
                                          nm_per_pixel * 0.001,
                                          wavelength * 0.001,
                                          refractive_index,
                                          numerical_aperture)

        self.pf = self.geo.createFromZernike(1.0, zmn)


class PupilFunctionScalar(PupilFunctionBase):
    """
    This version supports OTF scaling and sample aberrations.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, nm_per_pixel, zmn,
                 wavelength = pf_wavelength,
                 refractive_index = pf_refractive_index,
                 numerical_aperture = pf_numerical_aperture,
                 pf_size = pf_size,
                 otf_sigma = None,
                 sample_index = None,
                 z_center = None):
        """
        zmn is a list of lists containing the zernike mode terms, e.g.
            [[1.3, 2, 2]] for pure astigmatism.
        wavelength is the mean emission wavelength in nm.

        otf_sigma is the sigma to use for Gaussian OTF scaling.
        sample_index is the index of the sample media.
        z_center is the focal plane position in microns.

        If you specify sample_index you also need to specify z_center.
        """
        super(PupilFunctionScalar, self).__init__(sim_fp, x_size, y_size, h5_data, nm_per_pixel, pf_size)
        self.saveJSON({"psf" : {"class" : "PupilFunctionScalar",
                                "nm_per_pixel" : str(nm_per_pixel),
                                "numerical_aperture" : str(numerical_aperture),
                                "otf_sigma" : str(otf_sigma),
                                "pf_size" : str(pf_size),
                                "refactrive_index" : str(refractive_index),
                                "sample_index" : str(sample_index),
                                "wavelength" : str(wavelength),
                                "z_center" : str(z_center),
                                "zmn" : str(zmn)}})
        
        self.geo = pupilMath.Geometry(pf_size,
                                      nm_per_pixel * 0.001,
                                      wavelength * 0.001,
                                      refractive_index,
                                      numerical_aperture)

        self.pf = self.geo.createFromZernike(1.0, zmn)

        if otf_sigma is not None:
            self.otf_scaler = self.geo.gaussianScalingFactor(otf_sigma)
            
        self.sample_index = sample_index
        self.z_center = z_center

    def getAberrationFn(self, z_value):
        
        # Apply sample index aberration if requested.
        if self.sample_index is not None:

            # Distance 0.0 is the coverslip. Z values are relative to the focal
            # plane which is at z_center.
            #
            depth = self.z_center + z_value
            
            if (depth < 0.0):
                depth = 0.0

            return self.geo.aberration(depth, self.sample_index)

        else:
            return None
                

class PupilFunctionScalarCalibration(PupilFunctionBase):
    """
    This version supports OTF scaling and a fixed sample aberration. It
    simulates the PSF you'd expect from a calibration movie of a bead
    on a coverslip.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, nm_per_pixel, zmn,
                 wavelength = pf_wavelength,
                 refractive_index = pf_refractive_index,
                 numerical_aperture = pf_numerical_aperture,
                 pf_size = pf_size,
                 bead_z_center = None,
                 otf_sigma = None,
                 sample_index = None):
        """
        zmn is a list of lists containing the zernike mode terms, e.g.
            [[1.3, 2, 2]] for pure astigmatism.
        wavelength is the mean emission wavelength in nm.

        bead_z_center is the height of the bead above the coverslip in microns.
        otf_sigma is the sigma to use for Gaussian OTF scaling.
        sample_index is the index of the sample media.

        If you specify sample_index you also need to specify z_center.
        """
        super(PupilFunctionScalarCalibration, self).__init__(sim_fp, x_size, y_size, h5_data, nm_per_pixel, pf_size)
        self.saveJSON({"psf" : {"class" : "PupilFunctionScalarCalibration",
                                "bead_z_center" : str(bead_z_center),
                                "nm_per_pixel" : str(nm_per_pixel),
                                "numerical_aperture" : str(numerical_aperture),
                                "otf_sigma" : str(otf_sigma),
                                "pf_size" : str(pf_size),
                                "refactrive_index" : str(refractive_index),
                                "sample_index" : str(sample_index),
                                "wavelength" : str(wavelength),
                                "zmn" : str(zmn)}})
        
        self.geo = pupilMath.Geometry(pf_size,
                                      nm_per_pixel * 0.001,
                                      wavelength * 0.001,
                                      refractive_index,
                                      numerical_aperture)

        self.pf = self.geo.createFromZernike(1.0, zmn)

        if otf_sigma is not None:
            self.otf_scaler = self.geo.gaussianScalingFactor(otf_sigma)

        self.ab_fn = None
        if sample_index is not None:
            self.ab_fn = self.geo.aberration(bead_z_center, sample_index)

    def getAberrationFn(self, z_value):
        return self.ab_fn
        

class PupilFunctionVectorial(PupilFunctionBase):
    """
    This is a vectorial PSF version and supports bead scaling and sample aberrations.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, nm_per_pixel, zmn,
                 wavelength = pf_wavelength,
                 refractive_index = pf_refractive_index,
                 numerical_aperture = pf_numerical_aperture,
                 pf_size = pf_size,
                 bead_diameter = None,
                 sample_index = None,
                 z_center = None):
        """
        zmn is a list of lists containing the zernike mode terms, e.g.
            [[1.3, 2, 2]] for pure astigmatism.
        wavelength is the mean emission wavelength in nm.

        bead_diameter is the bead diameter in microns.
        sample_index is the index of the sample media.
        z_center is the focal plane position in microns.

        If you specify sample_index you also need to specify z_center.
        """
        super(PupilFunctionVectorial, self).__init__(sim_fp, x_size, y_size, h5_data, nm_per_pixel, pf_size)
        self.saveJSON({"psf" : {"class" : "PupilFunctionVectorial",
                                "nm_per_pixel" : str(nm_per_pixel),
                                "numerical_aperture" : str(numerical_aperture),
                                "bead_diameter" : str(bead_diameter),
                                "pf_size" : str(pf_size),
                                "refactrive_index" : str(refractive_index),
                                "sample_index" : str(sample_index),
                                "wavelength" : str(wavelength),
                                "z_center" : str(z_center),
                                "zmn" : str(zmn)}})
        
        self.geo = pupilMath.GeometryVectorial(pf_size,
                                               nm_per_pixel * 0.001,
                                               wavelength * 0.001,
                                               refractive_index,
                                               numerical_aperture)

        self.pf = self.geo.createFromZernike(1.0, zmn)

        if bead_diameter is not None:
            self.otf_scaler = self.geo.beadScalingFactor(bead_diameter)
            
        self.sample_index = sample_index
        self.z_center = z_center

    def getAberrationFn(self, z_value):
        
        # Apply sample index aberration if requested.
        if self.sample_index is not None:

            # Distance 0.0 is the coverslip. Z values are relative to the focal
            # plane which is at z_center.
            #
            depth = self.z_center + z_value
            
            if (depth < 0.0):
                depth = 0.0

            return self.geo.aberration(depth, self.sample_index)

        else:
            return None
        

class Spline2D(splineToPSF.SplineToPSF2D):
    """
    2D spline with non-zero offsets in x, y.
    """
    def __init__(self, spline_file):
        super(Spline2D, self).__init__(spline_file)
        self.psf_size = self.spline_size - 1

    def getPSF(self, z_value, dx, dy):
        """
        Use C library for better performance.

        z_value is ignored, it is only a parameter so that
        the signature matches Spline3D.

        dx, dy are in the range 0.0 - 1.0.
        """
        if((self.psf_size%2) == 0):
            return self.c_spline.getPSF(dy + 0.5, dx + 0.5)
        else:
            return self.c_spline.getPSF(dy + 1.0, dx + 1.0)

    def getPSFPy(self, z_value, dx, dy):
        """
        Python reference version.

        z_value is ignored, it is only a parameter so that
        the signature matches Spline3D.

        dx, dy are in the range 0.0 - 1.0.
        """
        psf = numpy.zeros((self.psf_size, self.psf_size))
        if((self.psf_size%2) == 0):
            for x in range(self.psf_size):
                for y in range(self.psf_size):
                    psf[y,x] = self.spline.f(float(y+dy) + 0.5,
                                             float(x+dx) + 0.5)
        else:
            for x in range(self.psf_size):
                for y in range(self.psf_size):
                    psf[y,x] = self.spline.f(float(y+dy) + 1.0,
                                             float(x+dx) + 1.0)
            
        return psf


class Spline3D(splineToPSF.SplineToPSF3D):
    """
    3D spline with non-zero offsets in x, y.
    """
    def __init__(self, spline_file):
        super(Spline3D, self).__init__(spline_file)
        self.psf_size = self.spline_size - 1

    def getPSF(self, z_value, dx, dy):
        """
        Use C library for better performance.

        z_value is ignored, it is only a parameter so that
        the signature matches Spline3D.

        dx, dy are in the range 0.0 - 1.0.
        """
        scaled_z = self.getScaledZ(z_value)
        if((self.psf_size%2) == 0):
            return self.c_spline.getPSF(scaled_z, dy + 0.5, dx + 0.5)
        else:
            return self.c_spline.getPSF(scaled_z, dy + 1.0, dx + 1.0)

    def getPSFPy(self, z_value, dx, dy):
        """
        Python reference version.

        z_value needs to be inside the z range covered by the spline.
        z_value should be in nanometers.

        dx, dy are in the range 0.0 - 1.0.
        """
        scaled_z = self.getScaledZ(z_value)

        psf = numpy.zeros((self.psf_size, self.psf_size))
        if((self.psf_size%2) == 0):
            for x in range(self.psf_size):
                for y in range(self.psf_size):
                    psf[y,x] = self.spline.f(scaled_z,
                                             float(y+dy) + 0.5,
                                             float(x+dx) + 0.5)
        else:
            for x in range(self.psf_size):
                for y in range(self.psf_size):
                    psf[y,x] = self.spline.f(scaled_z,
                                             float(y+dy) + 1.0,
                                             float(x+dx) + 1.0)

        return psf

    
class Spline(PSF):
    """
    PSF from a (cubic) spline.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, nm_per_pixel, spline_file):
        """
        spline_file is the name of a .spline file as generated by spliner/psf_to_spline.py.

        Splines are always pixel based, so the spline pixel size should match what you
        want for the simulation. The pixel size of the spline is recorded in the spline
        file and in theory that could be used to adjust appropriately depending on the
        desired pixel size for the simulation, but that is not currently being done.
        """
        super(Spline, self).__init__(sim_fp, x_size, y_size, h5_data, nm_per_pixel)
        self.saveJSON({"psf" : {"class" : "Spline",
                                "spline_file" : spline_file}})
        with open(spline_file, 'rb') as fp:
            spline_data = pickle.load(fp)
            
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

    def getPSFs(self, h5_data):
        """
        The expected form for the h5 data fields are x,y in pixels and z in microns.
        """
        image = numpy.zeros((self.im_size_x, self.im_size_y))
        x = h5_data['x']+1       # Pixels
        y = h5_data['y']+1       # Pixels
        z = 1.0*h5_data['z']*1000.0
        a = h5_data['sum']

        h5_data['height'] = numpy.zeros(a.size)

        dx = 1.0 - (x - numpy.floor(x))
        dy = 1.0 - (y - numpy.floor(y))

        for i in range(x.size):

            ix = int(x[i])
            iy = int(y[i])
            
            if (ix >= 0.0) and (ix < self.x_size) and (iy >= 0.0) and (iy < self.y_size):

                # Calculate psf, dx and dy are transposed to match the spline coordinate system.
                psf = self.spline.getPSF(z[i], dx[i], dy[i]).transpose()

                # Remove negative values.
                psf[(psf < 0.0)] = 0.0

                # Scale to correct number of photons.
                psf = a[i] * psf/numpy.sum(psf)
                h5_data['height'][i] = numpy.max(psf)

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
