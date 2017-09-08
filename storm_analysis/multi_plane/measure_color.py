#!/usr/bin/env python
"""
Performs all the steps necessary for color determination
from multi-color data.

Hazen 08/17
"""

import storm_analysis.multi_plane.copy_tracking as copyTracking
import storm_analysis.multi_plane.merge_heights as mergeHeights
import storm_analysis.multi_plane.mp_utilities_c as mpUtilC

import storm_analysis.sa_library.parameters as params

import storm_analysis.sa_utilities.avemlist_c as avemlistC


def measureColor(basename, n_planes):

    ref_name = basename + "mlist.bin"
    ch_alist_names = [basename + "alist.bin"]
    
    for i in range(1, n_planes):

        #
        # 1. Create localization files for the other channels with
        #    channel 0 tracking information.
        #
        ch_name = basename + "mlist_ch" + str(i) + ".bin"
        tracked_name = basename + "ch" + str(i) + "_tracked.bin"
        copyTracking.copyTracking(ref_name, ch_name, tracked_name)

        #
        # 2. Create averaged localization files for the tracked
        #    data for the other channels.
        #
        alist_name = basename + "alist_ch" + str(i) + ".bin"
        avemlistC.avemlist(tracked_name, alist_name)

        ch_alist_names.append(alist_name)

    #
    # 3. Merge height information from the different channels into
    #    a single file.
    #
    heights_name = basename + "heights.npy"
    mergeHeights.mergeHeights(heights_name, ch_alist_names)


if (__name__ == "__main__"):
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Measure heights from different channels.')

    parser.add_argument('--basename', dest='basename', type=str, required=True,
                        help = "The basename of the movie.")
    parser.add_argument('--xml', dest='xml', type=str, required=True,
                        help = "The analysis xml for this movie.")

    args = parser.parse_args()

    parameters = params.ParametersMultiplane().initFromFile(args.xml)

    measureColor(args.basename, mpUtilC.getNChannels(parameters))
