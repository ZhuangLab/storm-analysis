#!/usr/bin/env python
"""
Plot the channel height information from a data set.

Hazen 01/18
"""

import matplotlib
import matplotlib.pyplot as pyplot
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py

def plotHeights(h5_filename, channel_order, min_sum = 0.0):
    """
    h5_filename - The name of the HDF5 file.
    channel_order - A list of integers specifying the channel (wavelength) 
                    order, for example [0, 2, 1, 3]

    Note: The HDF5 file must have tracking information.
    """
    with saH5Py.SAH5Py(h5_filename) as h5:
        assert (h5.getNTracks() > 0), "No tracking information."
        assert (len(channel_order) == h5.getNChannels())

        # Figure out (ordered) height names.
        fields = []
        for elt in channel_order:
            if (elt == 0):
                fields.append("height")
            else:
                fields.append(h5.getChannelPrefix(elt) + "height")
                
        # Calculate total height.
        for tracks in h5.tracksIterator(fields = fields):
            total = numpy.zeros(tracks["height"].size)

            for i in range(h5.getNChannels()):
                total += tracks[fields[i]]

            # Only load the first block of tracks.
            break

        # Plot up to 100 points.
        n_pts = 0
        x = numpy.arange(h5.getNChannels())
        while(n_pts < 100):
            if (total[n_pts] > min_sum):
                y = []
                for i in range(h5.getNChannels()):
                    y.append(tracks[fields[i]][n_pts])
                y = numpy.array(y)
                pyplot.scatter(x + numpy.random.uniform(-0.2, 0.2),
                               y/total[n_pts],
                               color = "blue")
            n_pts += 1

        # Add labels.
        pyplot.ylim((0.0, 1.0))
        pyplot.xlabel("Channel Number")
        pyplot.ylabel("Normalized Height (AU)")
        pyplot.show()
    

if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Plots normalized channel heights.')

    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "Localization HDF5 file.")
    parser.add_argument('--order', nargs = '*', dest='order', type=str, required=True,
                        help = "The channel order by wavelength for example '2 1 0 3'.")

    args = parser.parse_args()
    
    channel_order = list(map(int, args.order))
    plotHeights(args.hdf5, channel_order)
