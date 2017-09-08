#!/usr/bin/env python
"""
Useful for reducing the imaging area and/or the number of
frames in a localization file.

Hazen 08/17
"""
import numpy
import sys
from xml.etree import ElementTree

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

def reduceMList(i3_name_in, i3_name_out, xstart = None, xstop = None, ystart = None, ystop = None, min_frame = None, max_frame = None):

    meta_data = readinsight3.loadI3Metadata(i3_name_in)
    
    i3_in = readinsight3.I3Reader(i3_name_in)
    i3_out = writeinsight3.I3Writer(i3_name_out)

    i3_data = i3_in.nextBlock()
    while (i3_data is not False):        
        sys.stdout.write(".")
        sys.stdout.flush()

        # Create mask.
        mask = numpy.full(i3_data.size, True, dtype = bool)
        if xstart is not None:
            mask = mask & (i3_data['xc'] > xstart)
        if xstop is not None:
            mask = mask & (i3_data['xc'] < xstop)
        if ystart is not None:
            mask = mask & (i3_data['yc'] > ystart)
        if ystop is not None:
            mask = mask & (i3_data['yc'] < ystop)
        if min_frame is not None:
            mask = mask & (i3_data['fr'] > min_frame)
        if max_frame is not None:
            mask = mask & (i3_data['fr'] < max_frame)

        i3_data = i3dtype.maskData(i3_data, mask)
        if (i3_data.size > 0):
            i3_out.addMolecules(i3_data)

        i3_data = i3_in.nextBlock()

    print()

    if meta_data is not None:
        i3_out.closeWithMetadata(ElementTree.tostring(meta_data, 'ISO-8859-1'))
    else:
        i3_out.close()        


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Reduce localizations list.')

    parser.add_argument('--in_bin', dest='in_bin', type=str, required=True,
                        help = "The name of the input localization file.")
    parser.add_argument('--out_bin', dest='out_bin', type=str, required=True,
                        help = "The name of the output localization file.")
    parser.add_argument('--xstart', dest='xstart', type=float, required=False,
                        help = "Minimum x position.")
    parser.add_argument('--xstop', dest='xstop', type=float, required=False,
                        help = "Maximum x position.")
    parser.add_argument('--ystart', dest='ystart', type=float, required=False,
                        help = "Minimum y position.")
    parser.add_argument('--ystop', dest='ystop', type=float, required=False,
                        help = "Maximum y position.")
    parser.add_argument('--min_frame', dest='min_frame', type=int, required=False,
                        help = "Minimum frame.")
    parser.add_argument('--max_frame', dest='max_frame', type=int, required=False,
                        help = "Maximum frame.")

    args = parser.parse_args()
        
    reduceMList(args.in_bin,
                args.out_bin,
                xstart = args.xstart,
                xstop = args.xstop,
                ystart = args.ystart,
                ystop = args.ystop,
                min_frame = args.min_frame,
                max_frame = args.max_frame)

