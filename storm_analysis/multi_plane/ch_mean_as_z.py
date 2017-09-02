#!/usr/bin/env python
"""
Given a list of paired localization files from a multi-color
experiment, create a localization file with mean channel
value as z and the total height in the 'a' field. Usually the
input files would be the alist files created by
multi_plane/batch_heights.py.

Why 'a'? Insight3 lets you filter localizations based on this
value.

Hazen 09/17
"""
import numpy

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

def ChannelMeanAsZ(input_filenames, output_filename):
    """
    Note: Input files should be ordered by channel wavelength.
    """

    # Create a reader for each file.
    i3_readers = []
    for i3_name in input_filenames:
        print(i3_name)
        i3_readers.append(readinsight3.I3Reader(i3_name))

    # Create writer for the results.
    i3_out = writeinsight3.I3Writer(output_filename)

    # Read first block of the first channel data.
    i3_data = [i3_readers[0].nextBlock()]
    while (i3_data[0] is not False):
        print("working..")

        # Read the data from the other channels.
        for i in range(1,len(i3_readers)):
            i3_data.append(i3_readers[i].nextBlock())

        # Calculate average and moment of channel heights.
        moment = numpy.zeros(i3_data[0].size)
        total = numpy.zeros(i3_data[0].size)

        for i in range(len(i3_readers)):
            h = i3_data[i]['h']
            print(" ", i, numpy.mean(h))
            moment += float(i)*h
            total += h

        moment = moment/total

        # Store moment and total in 'z' and 'a' fields respectively.
        i3_data[0]['z'] = 200.0*moment
        i3_data[0]['zc'] = 200.0*moment
        i3_data[0]['a'] = total

        i3_out.addMolecules(i3_data[0])

        # Load the next block of data.
        i3_data = [i3_readers[0].nextBlock()]

    # Close output file
    i3_out.close()

    
if (__name__ == "__main__"):
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Calculate color channel mean and store in localizations z value.')

    parser.add_argument('--alist', dest='alist', type=str, required=True, nargs = '*',
                        help = "Localization files in ordered by channel wavelength ")
    parser.add_argument('--output', dest='output', type=str, required=True,
                        help = "The name of the file for the results.")

    args = parser.parse_args()

    ChannelMeanAsZ(args.alist, args.output)
