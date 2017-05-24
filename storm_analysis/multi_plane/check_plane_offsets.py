#!/usr/bin/env python
"""
Given the measured PSF for each channel, plots their
maximum intensity as a function of z. This provides
some feedback on the actual spacing in z between the
different image planes.

Hazen 05/17
"""

import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import pickle


def checkPlaneOffsets(psf_files):

    psfs = []
    for psf_file in psf_files:
        with open(psf_file, 'rb') as fp:
            psfs.append(pickle.load(fp))
        
    fig = pyplot.figure()

    for psf in psfs:
        ave_psf = psf["psf"]
        zvals = psf["zvals"]

        max_i = numpy.max(ave_psf, axis = (1,2))
        pyplot.plot(zvals, max_i)

    pyplot.xlabel("Z offset (nm)")
    pyplot.ylabel("PSF Max (e-)")
    pyplot.show()


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Plots PSF max intensity versus for diagnostic purposes.')

    parser.add_argument('--psfs', nargs = '*', dest='psfs', type=str, required=True,
                        help = "The names of the PSF files to use.")

    args = parser.parse_args()

    checkPlaneOffsets(args.psfs)

    
