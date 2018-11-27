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

    plots = []
    for i, psf in enumerate(psfs):
        ave_psf = psf["psf"]
        zvals = psf["zvals"]

        max_i = numpy.amax(ave_psf, axis = (1,2))
        tmp, = pyplot.plot(zvals, max_i, label = psf_files[i])
        plots.append(tmp)

        n_max = numpy.argmax(max_i)
        print("Plane: {0:0d} maximum at {1:.1f}nm".format(i, zvals[n_max]))

    pyplot.legend(handles = plots, loc = 1)
    pyplot.ylim((0.0,1.5))
    pyplot.xlabel("Z offset (nm)")
    pyplot.ylabel("PSF Max (AU)")
    pyplot.show()


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Plots PSF max intensity versus for diagnostic purposes.')

    parser.add_argument('--psfs', nargs = '*', dest='psfs', type=str, required=True,
                        help = "The names of the PSF files to use.")

    args = parser.parse_args()

    checkPlaneOffsets(args.psfs)

    
