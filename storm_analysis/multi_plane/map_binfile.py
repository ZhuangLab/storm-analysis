#!/usr/bin/env python
"""
Map a i3 file to all channels. A debugging tool.

Hazen 05/17
"""

import os
import pickle

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3


def mapBinFile(i3_in_name, i3_out_name, xt, yt):
    i3_in = readinsight3.loadI3File(i3_in_name)
    
    xi = i3_in["x"]
    yi = i3_in["y"]

    xf = xt[0] + xt[1] * xi + xt[2] * yi
    yf = yt[0] + yt[1] * xi + yt[2] * yi
    
    i3dtype.posSet(i3_in, "x", xf)
    i3dtype.posSet(i3_in, "y", yf)

    with writeinsight3.I3Writer(i3_out_name) as wi:
        wi.addMolecules(i3_in)


if (__name__ == "__main__"):
    
    import argparse


    parser = argparse.ArgumentParser(description = 'Map a channel 0 i3 file to the other channels.')

    parser.add_argument('--bin', dest='i3bin', type=str, required=True,
                        help = "Channel0 Insight3 format binary file.")
    parser.add_argument('--mapping', dest='mapping', type=str, required=True,
                        help = "The name of the mapping file.")

    args = parser.parse_args()

    basename = os.path.splitext(args.i3bin)[0]

    # Load mappings.
    with open(args.mapping, 'rb') as fp:
        mappings = pickle.load(fp)

    for i in range(7):
        xt_name = "0_" + str(i+1) + "_x"
        yt_name = "0_" + str(i+1) + "_y"        
        if xt_name in mappings:
            xt = mappings[xt_name]
            yt = mappings[yt_name]
            mapBinFile(args.i3bin,
                       basename + "_c" + str(i+1) + ".bin",
                       xt,
                       yt)
