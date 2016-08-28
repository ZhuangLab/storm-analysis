#!/usr/bin/python
#
# Simple tagged spot file reader, to see if the file
# that I wrote is readable.
#
# Hazen 08/13
#

import sys

import google.protobuf.internal.decoder as decoder

import sa_library.readinsight3 as readinsight3
import TSFProto_pb2

if (len(sys.argv)!=2):
    print "usage: <tsf file>"
    exit()

tsf_file = open(sys.argv[1], "rb")

# Read magic number and offset
mnumber = readinsight3._getV(tsf_file, "I", 4)
offset = readinsight3._getV(tsf_file, ">Q", 8)

print "mnumber:", mnumber
print "offset:", offset

# Read SpotList message
tsf_file.seek(offset+12)

buffer = tsf_file.read()
(spot_list_size, position) = decoder._DecodeVarint(buffer, 0)

spot_list = TSFProto_pb2.SpotList()
spot_list.ParseFromString(buffer[position:position+spot_list_size])

print "Spot List:", spot_list_size
print spot_list
print ""

# Read the first few spots.
cur_pos = 12
for i in range(2):
    tsf_file.seek(cur_pos)

    buffer = tsf_file.read(200)
    (spot_size, position) = decoder._DecodeVarint(buffer, 0)

    spot = TSFProto_pb2.Spot()
    spot.ParseFromString(buffer[position:position+spot_size])

    cur_pos += position + spot_size

    print "Spot:", i
    print spot
    print ""


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
