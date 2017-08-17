#!/usr/bin/env python
"""
Map a i3 file to all channels, mostly a debugging tool.

Hazen 08/17
"""

import os
import pickle

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3


def mapXYLocations(xi, yi, xt, yt):
    xf = xt[0] + xt[1] * xi + xt[2] * yi
    yf = yt[0] + yt[1] * xi + yt[2] * yi
    return [xf, yf]
    

if (__name__ == "__main__"):
    
    import argparse


    parser = argparse.ArgumentParser(description = 'Map a channel0 i3 binary file to the other channels.')

    parser.add_argument('--bin', dest='i3bin', type=str, required=True,
                        help = "Localization file for a particular channel. Default is channel 0.")
    parser.add_argument('--mapping', dest='mapping', type=str, required=True,
                        help = "The name of the mapping file.")
    parser.add_argument('--channel', dest='channel', type=int, required=False, default=0,
                        help = "Channel that the localizations come from. Default is channel 0.")

    args = parser.parse_args()

    basename = os.path.splitext(args.i3bin)[0]

    # Load mappings.
    with open(args.mapping, 'rb') as fp:
        mappings = pickle.load(fp)

    # Load x,y locations.
    i3_in = readinsight3.loadI3File(args.i3bin)
    
    xi = i3_in["x"]
    yi = i3_in["y"]

    # Map to back channel 0 if necessary.
    if (args.channel != 0):
        xt_name = str(args.channel) + "_0_x"
        yt_name = str(args.channel) + "_0_y"
        [xi, yi] = mapXYLocations(xi, yi, mappings[xt_name], mappings[yt_name])

    for i in range(8):
        if (i == 0):
            i3dtype.posSet(i3_in, "x", xi)
            i3dtype.posSet(i3_in, "y", yi)
        else:
            xt_name = "0_" + str(i) + "_x"
            yt_name = "0_" + str(i) + "_y"        
            if xt_name in mappings:
                [xf, yf] = mapXYLocations(xi, yi, mappings[xt_name], mappings[yt_name])
                i3dtype.posSet(i3_in, "x", xf)
                i3dtype.posSet(i3_in, "y", yf)
            else:
                break

        i3_out_name = basename + "_mapb_c" + str(i) + ".bin"
        with writeinsight3.I3Writer(i3_out_name) as i3w:
            i3w.addMolecules(i3_in)

