#!/usr/bin/env python
"""
This is mostly for debugging. It takes the original molecule
list and creates one for each channel.

Note: This is not designed to be used with really large
      molecule list files.

Hazen 07/17
"""
import numpy
import os
import pickle

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

import storm_analysis.multi_plane.mp_utilities as mpUtil


def splitPeaks(mlist_filename, params_filename):

    parameters = params.ParametersMultiplane().initFromFile(params_filename)
    
    # Load the plane to plane mapping data.
    mappings = {}
    if parameters.hasAttr("mapping"):
        if os.path.exists(parameters.getAttr("mapping")):
            with open(parameters.getAttr("mapping"), 'rb') as fp:
                mappings = pickle.load(fp)

        else:
            print("No mapping file parameter, nothing to do")

    # Load frame offset information.
    frame_offsets = list(map(parameters.getAttr, mpUtil.getOffsetAttrs(parameters)))
    print(frame_offsets)
                  
    # Load the molecule list.
    ch0_data = readinsight3.loadI3File(mlist_filename)

    # Map to other channels.
    basename = mlist_filename[:-4]
    channel = 1
    m_key = "0_1_"
    while (m_key + "x") in mappings:
        chn_data = ch0_data.copy()

        # Map x.
        tx = mappings[m_key + "x"]
        chn_data['x'] = tx[0] + ch0_data['x']*tx[1] + ch0_data['y']*tx[2]

        # Map y.
        ty = mappings[m_key + "y"]
        chn_data['y'] = ty[0] + ch0_data['x']*ty[1] + ch0_data['y']*ty[2]

        # Map frame.
        chn_data['fr'] = ch0_data['fr'] - frame_offsets[0] + frame_offsets[channel]

        with writeinsight3.I3Writer(basename + "_ch" + str(channel) + ".bin") as i3w:
            i3w.addMolecules(chn_data)
                                    
        channel += 1
        m_key = "0_" + str(channel) + "_"
        

if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Split peaks')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file to split.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    splitPeaks(args.mlist, args.settings)
