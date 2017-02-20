#!/usr/bin/env python
"""
Checks that the on and off time distributions are as expected.

Hazen 12/16
"""

import argparse
import matplotlib
import matplotlib.pyplot as pyplot
import numpy

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3

parser = argparse.ArgumentParser(description = "Check emitter photophysics.")

parser.add_argument('--bin', dest='i3bin', type=str, required=True,
                    help = "The name of Insight3 format file with the emitter positions.")
parser.add_argument('--maxi', dest='maxi', type=int, required=False,
                    help = "The maximum number of emitters to check.")

args = parser.parse_args()

# Initialization.
i3_data = readinsight3.loadI3File(args.i3bin)
if args.maxi is None:
    args.maxi = numpy.max(i3_data['i']) + 1

if (args.maxi > (numpy.max(i3_data['i']) + 1)):
    args.maxi = numpy.max(i3_data['i']) + 1

max_frame = numpy.max(i3_data['fr'])

on_times = []
off_times = []
for i in range(int(args.maxi)):

    if ((i % 100) == 0):
        print("Analyzing emitter", i)

    # We are counting on getting the frames where the emitter was on in order.
    emitter_i3 = i3dtype.maskData(i3_data, (i3_data['i'] == i))
    frames_on = emitter_i3['fr']
    #print(frames_on, frames_on.size)

    # Check if the emitter was never on.
    if (frames_on.size == 0):
        continue

    # Check if off in the first frame.
    if (frames_on[0] != 1):
        off_times.append(frames_on[0])

    # Check if only a single frame.
    if (frames_on.size == 1):
        on_times.append(1)
        continue

    j = 1
    on_time = 1
    while (j < frames_on.size):
        cur_frame = frames_on[j]
        last_frame = frames_on[j-1]
        if (last_frame == (cur_frame - 1)):
            on_time += 1
        else:
            on_times.append(on_time)
            on_time = 1
            off_times.append(cur_frame - last_frame)
        j += 1
    if (cur_frame < max_frame):
        on_times.append(on_time)

if (len(on_times) > 0):
    on_times = numpy.array(on_times)
    print("On time mean", numpy.mean(on_times), "std", numpy.std(on_times))

    # histogram of on times.
    [counts, edges] = numpy.histogram(on_times, bins = numpy.arange(-0.5,10.5))
    centers = 0.5 * (edges[1:] + edges[:-1])
    fig = pyplot.figure()
    pyplot.plot(centers, counts)
    pyplot.xlabel("On Times (Frames)")
    pyplot.ylabel("Counts")
    pyplot.show()


if (len(off_times) > 0):
    off_times = numpy.array(off_times)
    print("Off time mean", numpy.mean(off_times), "std", numpy.std(off_times))

    # histogram of off times.
    [counts, edges] = numpy.histogram(off_times, bins = numpy.arange(-0.5,100.5,1))
    centers = 0.5 * (edges[1:] + edges[:-1])
    fig = pyplot.figure()
    pyplot.plot(centers, counts)
    pyplot.xlabel("Off Times (Frames)")
    pyplot.ylabel("Counts")
    pyplot.show()


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
