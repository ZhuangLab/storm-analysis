#!/usr/bin/env python
"""
Normalize the PSFs for each channel.

Hazen 05/17
"""

import pickle
import numpy
import os


def normalizePSFs(psf_files):

    # Load PSFs.
    psfs = []
    for psf_file in psf_files:
        with open(psf_file, 'rb') as fp:
            psfs.append(pickle.load(fp))

    # Determine relative maxima.
    psf_maxs = numpy.zeros(len(psfs))
    for i, psf in enumerate(psfs):
        psf_maxs[i] = numpy.amax(psf["psf"])

    psf_maxs = numpy.ones(len(psfs))/numpy.amax(psf_maxs)

    # Normalize PSFs. The brightest PSF will now have a maximum value
    # of 1.0, and other PSFs will have proportionally lower values.
    #
    # Also save the normalizations so that we can figure out how to
    # properly weight the different planes.
    #
    # FIXME: What is the correct way to weight? The relative total sum
    #        of the PSF? The height is close enough? More robust? It
    #        is probably okay as long as all the PSFs have basically
    #        the same shape..
    #
    for i, psf in enumerate(psfs):
        psf["psf"] = psf["psf"] * psf_maxs[i]
        psf["normalization"] = numpy.amax(psf["psf"])
        fname = os.path.splitext(psf_files[i])[0] + "_normed.psf"
        print(fname, psf["normalization"])
        with open(fname, 'wb') as fp:
            pickle.dump(psf, fp)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Normalize PSFs for multi-plane fitting.')

    parser.add_argument('--psfs', nargs = '*', dest='psfs', type=str, required=True,
                        help = "The names of the PSF files to normalize.")

    args = parser.parse_args()

    normalizePSFs(args.psfs)
    
