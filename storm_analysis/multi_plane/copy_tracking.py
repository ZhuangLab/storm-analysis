#!/usr/bin/env python
"""
This is part of the solution to the problem that our binary
format does not have enough channels for multi-color data.
 
You use it to copy the tracking information from channel 0
into the localization binary files for the other channels.
Once this has been done these binary files can be run through
sa_utilities/avem_list_c.py to generate a localization file
that has been averaged in the *same* way as the one that
was created for channel 0.

Hazen 08/17
"""
import numpy
import os

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

def copyTracking(ref_track_filename, input_filename, output_filename):

    i3_ref = readinsight3.I3Reader(ref_track_filename)
    i3_in = readinsight3.I3Reader(input_filename)

    # Verify that these two localization files are actually pairs.
    assert (i3_ref.molecules == i3_in.molecules), "Input files do not match."

#    # Verify that the output file doesn't already exist.
#    assert not os.path.exists(output_filename), "Output file already exists."
    
    i3_out = writeinsight3.I3Writer(output_filename)

    i3_ref_data = i3_ref.nextBlock()
    while (i3_ref_data is not False):        
        i3_in_data = i3_in.nextBlock()

        # Copy tracking information from the reference data file.
        for field in ["tl", "lk", "fi"]:
            i3_in_data[field] = i3_ref_data[field]

        # Save merge of reference (channel 0) and current channel.
        i3_out.addMolecules(i3_in_data)
        
        # Load the next block of reference data.
        i3_ref_data = i3_ref.nextBlock()

    # At least for now we are not bothering to copy the meta-data.
    i3_out.close()


if (__name__ == "__main__"):
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Copy tracking data from channel 0 into channel N.')

    parser.add_argument('--ch0', dest='ch0', type=str, required=True,
                        help = "The name of the channel 0 (tracked) localization file.")
    parser.add_argument('--chN', dest='chN', type=str, required=True,
                        help = "The name of the channel N localization file.")
    parser.add_argument('--output', dest='output', type=str, required=True,
                        help = "The name of the file for the results.")

    args = parser.parse_args()

    copyTracking(args.ch0, args.chN, args.output)
