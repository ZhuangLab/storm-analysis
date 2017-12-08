#!/usr/bin/env python
"""
For examining pixel level bias in localization positions.

Hazen 11/17
"""

import numpy

import storm_analysis.sa_library.readinsight3 as readinsight3


def pixelBias(filename, n_bins = 1000, normalized = True, i3_field = "x"):
    i3_data = readinsight3.loadI3File(filename)

    xv = numpy.fmod(i3_data[i3_field], 1.0)

    [xp_hist, x_bins] = numpy.histogram(xv, bins = n_bins, range = (0.0, 1.0), normed = normalized)
    x_centers = 0.5 * (x_bins[1:] + x_bins[:-1])

    return [x_centers, xp_hist]


if (__name__ == "__main__"):

    import argparse
    import matplotlib
    import matplotlib.pyplot as pyplot

    parser = argparse.ArgumentParser(description = 'Check fitter results for pixel level biases.')

    parser.add_argument('--locs', dest='locs', type=str, required=True,
                        help = "The name of the localizations binary file.")
    parser.add_argument('--n_bins', dest='n_bins', type=int, required=False, default = 1000,
                        help = "The number of bins in the histogram.")
    
    args = parser.parse_args()
        
    [x_centers, xp_hist] = pixelBias(args.locs, args.n_bins)
    
    pyplot.plot(x_centers, xp_hist)
    pyplot.show()

