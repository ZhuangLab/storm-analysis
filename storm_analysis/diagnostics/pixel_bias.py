#!/usr/bin/env python
"""
For examining pixel level bias in localization positions.

Hazen 08/18
"""

import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py


def pixelBias(h5_name, n_bins = 1000, normalized = True, field = "x"):

    sum_hist = None
    with saH5Py.SAH5Py(h5_name) as h5:
        for fnum, locs in h5.localizationsIterator(fields = [field]):
    
            xv = numpy.fmod(locs[field], 1.0)

            [xv_hist, xv_bins] = numpy.histogram(xv, bins = n_bins, range = (0.0, 1.0), normed = normalized)

            if sum_hist is None:
                x_centers = 0.5 * (xv_bins[1:] + xv_bins[:-1])
                sum_hist = xv_hist
            else:
                sum_hist += xv_hist

    return [x_centers, sum_hist]


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

    end_bin_sum = numpy.sum(xp_hist[:5] + xp_hist[-5:])
    print("Fraction of total localizations in the first and last 5 bins is {0:.3e}".format(end_bin_sum/numpy.sum(xp_hist)))
    print("Ideal fraction is {0:.3e}".format(10.0/args.n_bins))
    
    pyplot.plot(x_centers, xp_hist)
    pyplot.show()

