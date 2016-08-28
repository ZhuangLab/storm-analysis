#!/usr/bin/python
#
# Generate simulated data (2D).
#
# Hazen 01/16
#

import math
import numpy
import random
import sys

import storm_analysis.sa_library.daxwriter as daxwriter
import drawgaussians as dg
import storm_analysis.sa_library.writeinsight3 as writeinsight3

if (len(sys.argv) != 5):
    print "usage: <dax> <bin> <number frames> <number objects>"
    exit()

# Peak height and width.
sigma = 1.0
intensity = 500.0

# Image size.
x_size = 256
y_size = 256

dax_data = daxwriter.DaxWriter(sys.argv[1], x_size, y_size)
i3_data = writeinsight3.I3Writer(sys.argv[2])
num_frames = int(sys.argv[3])
num_objects = int(sys.argv[4])

for i in range(num_frames):
    print "Generating frame:", i

    # Generate the objects.
    objects = numpy.zeros((num_objects,5))
    for j in range(num_objects):
        if 1:
            x_off = 0.05 * float(x_size) + 0.9 * float(x_size) * random.random()
            y_off = 0.05 * float(y_size) + 0.9 * float(y_size) * random.random()
            objects[j,:] = [x_off,
                            y_off,
                            intensity,
                            sigma,
                            sigma]
            
        if 0:
            x_off = 10.0 + 10.0*(j%23) + 2.0*i
            y_off = 10.0 + 10.0*math.floor(float(j)/23.0)
            objects[j,:] = [x_off,
                            y_off,
                            intensity,
                            sigma,
                            sigma]

    # Draw the image.
    #image = dg.drawGaussians([x_size, y_size], objects, background = 200, res = 5)
    image = dg.drawGaussians([x_size, y_size], objects, background = 0, res = 5)

    # Add poisson noise and baseline.
    image = numpy.random.poisson(image) + 100.0
    
    # Save the image
    dax_data.addFrame(image)

    # Save the molecule locations.
    frame = (i+1)*numpy.ones(objects.shape[0])
    i3_data.addMoleculesWithXYIFrame(objects[:,0] + 1.0,
                                     objects[:,1] + 1.0,
                                     objects[:,2],
                                     frame)

dax_data.close()
i3_data.close()


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
