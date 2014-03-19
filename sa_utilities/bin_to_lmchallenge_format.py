#!/usr/bin/python
#
# Convert a .bin format file to the Single Molecule
# Localization Challenge format.
#
# Hazen 03/14
#

import sys

import sa_library.readinsight3 as readinsight3

if (len(sys.argv)!=4):
    print "usage: <bin_file> <smlc_file> <pix_to_nm>"
    exit()

i3_reader = readinsight3.I3Reader(sys.argv[1])
i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

smlc_file_fp = open(sys.argv[2], "w")
smlc_file_fp.write("frame,xnano,ynano,intensity\n")

pix_to_nm = float(sys.argv[3])

print "Saving Localizations"
localization_number = 0
while (type(i3_block) != type(False)):

    print " saving localization", localization_number

    for i in range(len(i3_block)):
        fr = i3_block['fr'][i]
        yp = (i3_block['xc'][i] - 0.5) * pix_to_nm
        xp = (i3_block['yc'][i] - 0.5) * pix_to_nm
        intensity = i3_block['a'][i]

        smlc_file_fp.write("{0:d}, {1:.3f}, {2:.3f}, {3:.3f}\n".format(fr, xp, yp, intensity))

    localization_number += len(i3_block)
    i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

smlc_file_fp.close()
