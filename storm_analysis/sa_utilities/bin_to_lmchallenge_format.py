#!/usr/bin/env python
"""
03/14

Convert a .bin format file to the Single Molecule
Localization Challenge format.

07/16

Updated for the 2016 challenge.

Hazen
"""

import sys

import storm_analysis.sa_library.readinsight3 as readinsight3


def binToLMChallenge(bin_name, smlc_name, pix_to_nm, verbose = True):
    """
    A molecule that was on for several frames is stored as a single track,
    so each one is expanded back into one row per frame.

    'a' is the total number of photons over the whole track, both in the
    files this package writes and in the Insight3 era ones, where avemlist
    flagged AREA as a total. So each row gets its share of it rather than
    all of it.

    Note that a track with a gap in it is expanded into consecutive frames
    anyway, since the .bin format records only where the track started and
    how long it was.
    """
    i3_reader = readinsight3.I3Reader(bin_name)
    i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

    smlc_file_fp = open(smlc_name, "w")
    smlc_file_fp.write("index, frame, xnano, ynano, znano, intensity\n")

    if verbose:
        print("Saving Localizations")
    localization_number = 0
    index = 0
    while (type(i3_block) != type(False)):

        if verbose:
            print(" saving localization", localization_number)

        for i in range(len(i3_block)):
            track_length = int(i3_block['tl'][i])
            if (track_length < 1):
                track_length = 1
            for j in range(track_length):
                fr = i3_block['fr'][i] + j
                xp = i3_block['xc'][i] * pix_to_nm
                yp = i3_block['yc'][i] * pix_to_nm
                zp = i3_block['zc'][i]
                intensity = i3_block['a'][i] / float(track_length)

                index += 1
                smlc_file_fp.write("{0:d}, {1:d}, {2:.3f}, {3:.3f}, {4:.3f}, {5:.3f}\n".format(index, fr, xp, yp, zp, intensity))

        localization_number += len(i3_block)
        i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

    if verbose:
        print("Saved", index, "molecules.")
    smlc_file_fp.close()

    return index


if (__name__ == "__main__"):

    if (len(sys.argv)!=4):
        print("usage: <bin_file> <smlc_file> <pix_to_nm>")
        exit()

    binToLMChallenge(sys.argv[1], sys.argv[2], float(sys.argv[3]))

#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
