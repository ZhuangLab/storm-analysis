#!/usr/bin/env python
"""
Some math for calculating PSFs from pupil functions.

All units are in microns.

Hazen 03/16
"""

import math
import numpy
import scipy
import scipy.fftpack
import tifffile

import storm_analysis.simulator.zernike_c as zernikeC

#
# FIXME: By eye I would say that this is off by about a factor of two.
#        The PSF that this generates is about 2x narrower than expected.
#
#        I think is fixed now, or at least it looks more like what I
#        would expect. Still needs testing and some math skills..
#        (2016-12-07 HB).
#
class Geometry(object):

    def __init__(self, size, pixel_size, wavelength, imm_index, NA):
        """
        size - The number of pixels in the PSF image, assumed square.
        pixel_size - The size of the camera pixel in um.
        wavelength - The wavelength of the flourescence in um.
        imm_index - The index of the immersion media.
        NA - The numerical aperature of the objective.
        """

        self.imm_index = float(imm_index)
        self.NA = float(NA)
        self.pixel_size = float(pixel_size)
        self.size = int(size)
        self.wavelength = float(wavelength)

        self.k_max = NA/wavelength

        dk = 2.0/(size * pixel_size)
        self.r_max = self.k_max/dk
        
        [x,y] = numpy.mgrid[ -self.size/2.0 : self.size/2.0, -self.size/2.0 : self.size/2.0]
        kx = dk * x
        ky = dk * y
        self.k = numpy.sqrt(kx * kx + ky * ky)
        self.kx = kx * 2.0/size
        self.ky = ky * 2.0/size
        
        tmp = imm_index/wavelength
        self.kz = numpy.lib.scimath.sqrt(tmp * tmp - self.k * self.k)

        self.r = self.k/self.k_max

        self.kz[(self.r > 1.0)] = 0.0
        self.n_pixels = numpy.sum(self.r <= 1)
        self.norm = math.sqrt(self.r.size)

    def applyNARestriction(self, pupil_fn):
        """
        pupil_fn - The pupil function to restrict the NA of.
    
        return - The NA restricted pupil function.
        """
        pupil_fn[(self.r > 1.0)] = 0.0
        return pupil_fn

    def changeFocus(self, pupil_fn, z_dist):
        """
        pupil_fn - The pupil function.
        z_dist - The distance to the new focal plane.
        
        return - The pupil function at the new focal plane.
        """
        return numpy.exp(1j * 2.0 * numpy.pi * self.kz * z_dist) * pupil_fn

    def createPlaneWave(self, n_photons):
        """
        n_photons - The intensity of the pupil function.
    
        return - The pupil function for a plane wave.
        """
        plane = numpy.sqrt(n_photons/self.n_pixels) * numpy.exp(1j * numpy.zeros(self.r.shape))
        return self.applyNARestriction(plane)

    def createFromZernike(self, n_photons, zernike_modes):
        """
        n_photons - The intensity of the pupil function
        zernike_modes - List of lists, [[magnitude (in radians), m, n], [..]]
    
        return - The pupil function for this combination of zernike modes.
        """
        if (len(zernike_modes) == 0):
            return self.createPlaneWave(n_photons)
        else:
            phases = numpy.zeros(self.r.shape)
            for zmn in zernike_modes:
                phases = zernikeC.zernikeGrid(phases, zmn[0], zmn[1], zmn[2], radius = self.r_max)
            zmnpf = numpy.sqrt(n_photons/self.n_pixels) * numpy.exp(1j * phases)
            return self.applyNARestriction(zmnpf)

    def dx(self, pupil_fn):
        return -1j * 2.0 * numpy.pi * self.kx * pupil_fn
    
    def pfToPSF(self, pf, z_vals, want_intensity = True, scaling_factor = None):
        """
        pf - A pupil function.
        z_vals - The z values (focal planes) of the desired PSF.
        want_intensity - (Optional) Return intensity, default is True.
        scaling_factor - (Optional) The OTF rescaling factor, default is None.

        return - The PSF that corresponds to pf at the requested z_vals.
        """
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

    def translatePf(self, pupil_fn, dx, dy):
        """
        Translate the Pf using Fourier translation.
    
        pupil_fn - A pupil function.
        dx - Translation in x in pixels.
        dy - Translation in y in pixels.
    
        return - The PF translated by dx, dy.
        """
        return numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy)) * pupil_fn

    
def intensity(x):
    """
    x - The (numpy array) to convert to intensity.

    return - The product of x and the complex conjugate of x.
    """
    return numpy.abs(x * numpy.conj(x))


def toRealSpace(pupil_fn):
    """
    pupil_fn - A pupil function.

    return - The pupil function in real space (as opposed to fourier space).
    """
    return scipy.fftpack.ifftshift(math.sqrt(pupil_fn.size) * scipy.fftpack.ifft2(pupil_fn))


if (__name__ == "__main__"):

    import pickle
    import sys

    if (len(sys.argv) < 2):
        print("usage: <psf> <zmn.txt> <amp>")
        exit()

    pixel_size = 0.10
    #pixel_size = 0.020
    wavelength = 0.6
    refractive_index = 1.5
    numerical_aperture = 1.4
    z_range = 1.0
    z_pixel_size = 0.010

    geo = Geometry(int(20.0/pixel_size),
                   pixel_size,
                   wavelength,
                   refractive_index,
                   numerical_aperture)

    if (len(sys.argv) == 4):
        zmn = []
        amp = float(sys.argv[3])
        with open(sys.argv[2]) as fp:
            for line in fp:
                data = line.strip().split(" ")
                if (len(data) == 3):
                    zmn.append([amp * float(data[2]), int(data[0]), int(data[1])])
    else:
        zmn = [[1.3, 2, 2]]
        #zmn = [[1, 0, 4]]
        #zmn = []

    pf = geo.createFromZernike(1.0, zmn)
    z_values = numpy.arange(-z_range, z_range + 0.5 * z_pixel_size, z_pixel_size)
    psfs = geo.pfToPSF(pf, z_values)

    #xy_size = 2.0*psfs.shape[0]
    xy_size = 100
    xy_start = int(0.5 * (psfs.shape[1] - xy_size) + 1)
    xy_end = int(xy_start + xy_size)
    psfs = psfs[:,xy_start:xy_end,xy_start:xy_end]
    
    if 1:
        tifffile.imsave("pf_abs.tif", numpy.abs(pf).astype(numpy.float32))
        tifffile.imsave("pf_angle.tif", (180.0 * numpy.angle(pf)/numpy.pi + 180).astype(numpy.float32))

    if 1:
        with tifffile.TiffWriter(sys.argv[1]) as psf_tif:
            temp = (psfs/numpy.max(psfs)).astype(numpy.float32)
            psf_tif.save(temp)

    if 1:
        with open("z_offset.txt", "w") as fp:
            for i in range(z_values.size):
                fp.write("1 {0:.6f}\n".format(1000.0 * z_values[i]))
        
    if 0:
        psfs = (65000.0 * (psfs/numpy.max(psfs))).astype(numpy.uint16)
        psf_dict = {"pixel_size" : pixel_size,
                    "wavelength" : wavelength,
                    "refractive_index" : refractive_index,
                    "numerical_aperture" : numerical_aperture,
                    "z_range" : z_range,
                    "zmm" : zmn,
                    "psf" : psfs}
        pickle.dump(psf_dict, open(sys.argv[1], "wb"), protocol = 2)


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


