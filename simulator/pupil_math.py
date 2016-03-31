#!/usr/bin/python
#
# Some math for calculating PSFs from pupil functions.
#
# All units are in microns.
#
# Hazen 03/16
#

import math
import numpy
import scipy
import scipy.fftpack
import scipy.optimize
import sys
import tifffile

import zernike_c as zernikeC

class Geometry(object):

    ## __init__
    #
    # @param size The number of pixels in the PSF image, assumed square.
    # @param pixel_size The size of the camera pixel in nm.
    # @param wavelength The wavelength of the flourescence in um.
    # @param imm_index The index of the immersion media.
    # @param The numerical aperature of the objective.
    #
    def __init__(self, size, pixel_size, wavelength, imm_index, NA):

        self.imm_index = float(imm_index)
        self.NA = float(NA)
        self.pixel_size = float(pixel_size)
        self.size = int(size)
        self.wavelength = float(wavelength)

        self.k_max = NA/wavelength

        dk = 1.0/(size * pixel_size)
        self.r_max = self.k_max/dk
        
        [x,y] = numpy.mgrid[ -size/2.0 : size/2.0, -size/2.0 : size/2.0] + 0.5
        kx = dk * x
        ky = dk * y
        self.k = numpy.sqrt(kx*kx + ky*ky)
        
        tmp = imm_index/wavelength
        self.kz = numpy.lib.scimath.sqrt(tmp * tmp - self.k * self.k)

        self.r = self.k/self.k_max
        self.n_pixels = numpy.sum(self.r <= 1)
        self.norm = math.sqrt(self.r.size)

    ## applyNARestriction
    #
    # @param pupil_fn The pupil function to restrict the NA of.
    #
    # @return The NA restricted pupil function.
    #
    def applyNARestriction(self, pupil_fn):
        pupil_fn[(self.r > 1.0)] = 0.0
        return pupil_fn

    ## changeFocus
    #
    # @param pupil_fn The pupil function.
    # @param z_dist The distance to the new focal plane.
    #
    # @return The pupil function at the new focal plane.
    #
    def changeFocus(self, pupil_fn, z_dist):
        return numpy.exp(2.0 * numpy.pi * 1j * self.kz * z_dist) * pupil_fn

    ## createPlaneWave
    #
    # @param n_photons The intensity of the pupil function.
    #
    # @return The pupil function for a plane wave.
    #
    def createPlaneWave(self, n_photons):
        plane = numpy.sqrt(n_photons/self.n_pixels) * numpy.exp(1j * numpy.zeros(self.r.shape))
        return self.applyNARestriction(plane)

    ## createFromZernike
    #
    # @param n_photons The intensity of the pupil function
    # @param zernike_modes List of lists, [[magnitude (in radians), m, n], [..]]
    #
    # @return The pupil function for this combination of zernike modes.
    #
    def createFromZernike(self, n_photons, zernike_modes):
        phases = numpy.zeros(self.r.shape)
        for zmn in zernike_modes:
            phases = zernikeC.zernikeGrid(phases, zmn[0], zmn[1], zmn[2], radius = self.r_max)
        zmnpf = numpy.sqrt(n_photons/self.n_pixels) * numpy.exp(1j * phases)
        return self.applyNARestriction(zmnpf)
    
    ## pfToPSF
    #
    # @param pf A pupil function.
    # @param z_vals The z values (focal planes) of the desired PSF.
    # @param want_intensity (Optional) Return intensity, default is True.
    # @param scaling_factor (Optional) The OTF rescaling factor, default is None.
    #
    # @return The PSF that corresponds to pf at the requested z_vals.
    #
    def pfToPSF(self, pf, z_vals, want_intensity = True, scaling_factor = None):
        if want_intensity:
            psf = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]))
            for i, z in enumerate(z_vals):
                defocused = toRealSpace(self.changeFocus(pf, z))
                if scaling_factor is not None:
                    otf = scipy.fftpack.fftshift(scipy.fftpack.fft2(intensity(defocused)))
                    otf_scaled = otf * scaling_factor
                    psf[i,:,:] = numpy.abs(scipy.fftpack.ifft2(otf_scaled))
                else:
                    psf[i,:,:] = intensity(defocused)
            return psf
        else:
            psf = numpy.zeros((len(z_vals), pf.shape[0], pf.shape[1]),
                              dtype = numpy.complex_)
            for i, z in enumerate(z_vals):
                psf[i,:,:] = toRealSpace(self.changeFocus(pf, z))
            return psf


## intensity
#
# @param x The (numpy array) to convert to intensity.
#
# @return The product of x and the complex conjugate of x.
#
def intensity(x):
    return numpy.abs(x * numpy.conj(x))

## toRealSpace
#
# @param pupil_fn A pupil function.
#
# @return The pupil function in real space (as opposed to fourier space).
#
def toRealSpace(pupil_fn):
    return scipy.fftpack.ifftshift(math.sqrt(pupil_fn.size) * scipy.fftpack.ifft2(pupil_fn))


if (__name__ == "__main__"):
    
    geo = Geometry(1600, 0.004, 0.6, 1.5, 1.4)

    zmn = [[2.0, 2, 2]]
    #pf = geo.createPlaneWave(1.0)
    pf = geo.createFromZernike(1.0, zmn)

    if 1:
        tifffile.imsave("pf_abs.tif", (1000.0 * numpy.abs(pf)).astype(numpy.uint16))
        tifffile.imsave("pf_angle.tif", (1800.0 * numpy.angle(pf)/numpy.pi + 1800).astype(numpy.uint16))

    if 1:
        psfs = geo.pfToPSF(pf, [-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4])
        
        with tifffile.TiffWriter("psf.tif") as psf_tif:
            for i in range(psfs.shape[0]):
                temp = (psfs[i,:,:] * 1.0e6).astype(numpy.uint16)
                psf_tif.save(temp)


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


