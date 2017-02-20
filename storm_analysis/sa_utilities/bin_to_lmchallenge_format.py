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

if (len(sys.argv)!=4):
    print("usage: <bin_file> <smlc_file> <pix_to_nm>")
    exit()

i3_reader = readinsight3.I3Reader(sys.argv[1])
i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

smlc_file_fp = open(sys.argv[2], "w")
smlc_file_fp.write("index, frame, xnano, ynano, znano, intensity\n")

pix_to_nm = float(sys.argv[3])

print("Saving Localizations")
localization_number = 0
index = 0
while (type(i3_block) != type(False)):

    print(" saving localization", localization_number)

    for i in range(len(i3_block)):
        track_length = i3_block['tl'][i]
        for j in range(track_length):
            fr = i3_block['fr'][i] + j
            xp = i3_block['xc'][i] * pix_to_nm
            yp = i3_block['yc'][i] * pix_to_nm
            zp = i3_block['zc'][i]
            intensity = i3_block['a'][i]

            index += 1
            smlc_file_fp.write("{0:d}, {1:d}, {2:.3f}, {3:.3f}, {4:.3f}, {5:.3f}\n".format(index, fr, xp, yp, zp, intensity))

    localization_number += len(i3_block)
    i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

print("Saved", index, "molecules.")
smlc_file_fp.close()

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
