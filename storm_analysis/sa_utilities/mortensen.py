#!/usr/bin/env python
"""
Calculate CRao bound for localization using eq. 5 in 
Mortensen - Nature Methods - 2010.

Hazen 12/17
"""
import math
import numpy
import scipy
import scipy.integrate


def cramerRaoBound(intensity, background, pixel_size, psf_size, is_emccd = False):
    """
    intensity - photo-electrons.
    background - photo-electrons.
    pixel_size - camera pixel size in nanometers.
    psf_size - PSF sigma in nanometers.
    """
    
    px_sqr = pixel_size * pixel_size

    #
    # This is the average value returned by daostorm analysis, 2d
    # fit for the highest intensity bead data.
    #
    psf_sqr = psf_size * psf_size

    sa_sqr = psf_sqr + px_sqr/12.0

    def integral_fn(t, N, bg_sqr):
        ln_t = math.log(t)
        t2 = N * px_sqr * t / (2.0 * math.pi * sa_sqr * bg_sqr)
        return ln_t/(1.0 + t2)

    integ = scipy.integrate.quad(lambda x: integral_fn(x, intensity, background), 1.0e-9, 1.0)[0]

    if is_emccd:
        return math.sqrt(2.0 * (sa_sqr/intensity) * 1.0/(1.0 + integ))
    else:
        return math.sqrt((sa_sqr/intensity) * 1.0/(1.0 + integ))


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Cramer-Rao Bounds - Mortensen, Nature Methods, 2010')

    parser.add_argument('--intensity', dest='intensity', type=float, required=True,
                        help = "Localization intensity in photons.")
    parser.add_argument('--background', dest='background', type=float, required=True,
                        help = "Per pixel background in photons.")
    parser.add_argument('--pixel_size', dest='pixel_size', type=float, required=True,
                        help = "Pixel size in nanometers.")
    parser.add_argument('--psf_sigma', dest='psf_sigma', type=float, required=True,
                        help = "PSF sigma in nanometers.")

    args = parser.parse_args()
    
    print("X/Y bound {0:.3f} nm".format(cramerRaoBound(args.intensity,
                                                       args.background,
                                                       args.pixel_size,
                                                       args.psf_sigma)))
    
