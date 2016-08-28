#!/usr/bin/python
#
# Convert a .bin format file to the micro-manager tagged
# spot format file.
#
# Hazen 08/13
#

import os
import struct
import sys

import google.protobuf.internal.encoder as encoder

import sa_library.datareader as datareader
import sa_library.readinsight3 as readinsight3
import TSFProto_pb2

if (len(sys.argv)!=5):
    print "usage: <dax file> <bin file> <tsf file> <pixel size (nm)>"
    exit()

# Open TSF file.
#if os.path.exists(sys.argv[3]):
#    print "attempting to overwrite an existing tsf file, exitting"
#    exit()

def setV(fp, format, value):
    fp.write(struct.pack(format, value))

pix_to_nm = float(sys.argv[4])

tsf_file = open(sys.argv[3], "wb")
setV(tsf_file, "I", 0)
setV(tsf_file, ">Q", 0)

# Save Spot(s).
i3_reader = readinsight3.I3Reader(sys.argv[2])
i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

print "Saving localizations"
channels = []
localization_number = 0
while (type(i3_block) != type(False)):

    print " saving localization", localization_number

    for i in range(len(i3_block)):
        spot = TSFProto_pb2.Spot()

        spot.molecule = localization_number
        channel = int(i3_block['c'][i])
        if not channel in channels:
            channels.append(channel)
        spot.channel = channel
        spot.frame = int(i3_block['fr'][i])

        # These are always the same.
        spot.slice = 1
        spot.pos = 1

        spot.x = float(i3_block['xc'][i])*pix_to_nm
        spot.y = float(i3_block['yc'][i])*pix_to_nm
        spot.z = float(i3_block['zc'][i])

        spot.intensity = float(i3_block['a'][i])
        spot.background = float(i3_block['bg'][i])
        spot.width = float(i3_block['w'][i])
        spot.a = float(i3_block['ax'][i])
        spot.theta = 0.0

        spot.x_original = float(i3_block['x'][i])*pix_to_nm
        spot.y_original = float(i3_block['y'][i])*pix_to_nm
        spot.z_original = float(i3_block['z'][i])

        localization_number += 1

        out = spot.SerializeToString()
        out = encoder._VarintBytes(len(out)) + out
        tsf_file.write(out)

    i3_block = i3_reader.nextBlock(block_size = 1000, good_only = False)

print ""
print localization_number, "total localizations"
# Save SpotList.
print ""
print "Saving analysis meta-data"
print " data set contains", len(channels), "channels", channels
spot_list = TSFProto_pb2.SpotList()

# FIXME: get a real id..
spot_list.application_id = 1
spot_list.name = os.path.basename(sys.argv[1])
spot_list.filepath = spot_list.name
spot_list.pixel_size = pix_to_nm
spot_list.nr_spots = localization_number

spot_list.nr_channels = len(channels)

# These are always the same.
spot_list.nr_slices = 1
spot_list.nr_pos = 1
spot_list.fit_mode = 1
spot_list.location_units = 0
spot_list.intensity_units = 0
spot_list.is_track = False

# If a dax file is provided, get the film size.
if os.path.exists(sys.argv[1]):
    data_reader = datareader.inferReader(sys.argv[1])
    [x, y, l] = data_reader.filmSize()

    spot_list.nr_pixels_x = x
    spot_list.nr_pixels_y = y
    spot_list.nr_frames = l

spot_list_offset = tsf_file.tell() - 12

out = spot_list.SerializeToString()
out = encoder._VarintBytes(len(out)) + out
tsf_file.write(out)

# Rewind to the beginning and record the offset of the SpotList message.
tsf_file.seek(4)
setV(tsf_file, ">Q", spot_list_offset)
tsf_file.close()

#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
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
