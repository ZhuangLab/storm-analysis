#!/usr/bin/env python
"""
Plot the channel height information from a data set.

Hazen 08/17
"""

import matplotlib
import matplotlib.pyplot as pyplot
import numpy


def plotHeights(ch_order, heights, min_sum = 0.0):

    # Calculate height totals.
    total_h = numpy.sum(heights, axis = 1)

    # Remove localizations with small total height.
    if (min_sum > 0.0):
        mask = (total_h > min_sum)

        total_h = total_h[mask]
        heights = heights[mask,:]
        
    x = numpy.array(ch_order)
    y = heights

    print(x.size, y.shape)

    fig = pyplot.figure()

    # Plot up to 100 points.
    for i in range(100):
        if (i < y.shape[0]):
            pyplot.scatter(x + numpy.random.uniform(-0.2, 0.2), y[i,:]/total_h[i])

    # Add labels.
    pyplot.ylim((0.0, 1.0))
    pyplot.xlabel("Channel Number")
    pyplot.ylabel("Normalized Height (AU)")
    pyplot.show()
    

if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Plots normalized channel heights for diagnostic purposes.')

    parser.add_argument('--heights', dest='heights', type=str, required=True,
                        help = "The name of the heights file to plot.")
    parser.add_argument('--order', dest='order', type=int, required=True, nargs = '*',
                        help = "The channel order.")
    parser.add_argument('--min_sum', dest='min_sum', type=float, required=False, default = 0.0,
                        help = "The minimum heights sum. Default 0.0")

    args = parser.parse_args()

    height_data = numpy.load(args.heights)
    plotHeights(args.order, height_data, min_sum = args.min_sum)
