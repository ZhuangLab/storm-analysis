#!/usr/bin/env python
"""
This is part of the solution to the problem that our binary
format does not have enough channels for multi-color data.

You use it to merge the heights from each channel into a
single .npy file that is used by downstream analysis to
assign a color to the localization.

Hazen 08/17
"""
import numpy

import storm_analysis.sa_library.readinsight3 as readinsight3

def mergeHeights(height_filename, channel_filenames):
    """
    channel_filenames is a list in order from shortest wavelength
    to the longest wavelength.
    """
    #
    # FIXME: We are just loading all the data at once, which could
    #        problematic for really large files.
    #
    height_data = None
    for i, elt in enumerate(channel_filenames):
        i3_data = readinsight3.loadI3File(elt)

        if height_data is None:
            height_data = numpy.zeros((i3_data.size, len(channel_filenames)), dtype = numpy.float32)

        height_data[:,i] = i3_data['h']

    # User feedback.
    if (height_data.shape[0] > 5):
        for i in range(5):
            print(height_data[i,:])

    # Save the heights.
    numpy.save(height_filename, height_data)


if (__name__ == "__main__"):
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Merge heights from different channels.')

    parser.add_argument('--output', dest='output', type=str, required=True,
                        help = "The name of the (numpy) file to save the height data in.")
    parser.add_argument('--channels', dest='channels', type=str, required=True, nargs = '*',
                        help = "The names of the channel localization files in order from shortest to longest wavelength.")

    args = parser.parse_args()

    mergeHeights(args.output, args.channels)
    
        
