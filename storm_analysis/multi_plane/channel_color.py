#!/usr/bin/env python
"""
Given a tracked multiplane HDF5 file this will calculate the
total height for each track across all the channels as well
as the first moment of the height with respect to channel
number.

Hazen 01/18
"""
import copy
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py

def channelColor(h5_filename, channel_order):
    """
    h5_filename - The name of the HDF5 file.
    channel_order - A list of integers specifying the channel (wavelength) 
                    order, for example [0, 2, 1, 3]

    Note: The HDF5 file must have tracking information.
    """
    with saH5Py.SAH5Py(h5_filename) as h5:
        assert (h5.getNTracks() > 0), "No tracking information."
        assert (len(channel_order) == h5.getNChannels())

        # Figure out order height names.
        height_names = []
        for elt in channel_order:
            if (elt == 0):
                height_names.append("height")
            else:
                height_names.append(h5.getChannelPrefix(elt) + "height")
        
        fields = copy.copy(height_names)
            
        for i, tracks in enumerate(h5.tracksIterator(fields = fields)):
            moment = numpy.zeros(tracks["height"].size)
            total = numpy.zeros(tracks["height"].size)

            for j in range(h5.getNChannels()):
                moment += float(j) * tracks[height_names[j]]
                total += tracks[height_names[j]]

            moment = moment/total

            h5.addTrackData(moment, i, "height_moment")
            h5.addTrackData(total, i, "height_total")

    
if (__name__ == "__main__"):
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Calculate color channel mean and total.')

    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "Localization HDF5 file.")
    parser.add_argument('--order', nargs = '*', dest='order', type=str, required=True,
                        help = "The channel order by wavelength for example '2 1 0 3'.")

    args = parser.parse_args()
    
    channel_order = list(map(int, args.order))
    channelColor(args.hdf5, channel_order)
