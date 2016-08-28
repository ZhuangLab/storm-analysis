#!/usr/bin/python
#
# Generate simulated data (3D).
#
# Hazen 01/16
#

import math
import numpy
import random
import sys

import sa_library.daxwriter as daxwriter
import drawgaussians as dg
import sa_library.i3dtype as i3dtype
import sa_library.writeinsight3 as writeinsight3

#import astigmaticPSF as PSF
import dhPSF as PSF

if (len(sys.argv) != 5):
    print "usage: <dax> <bin> <frames> <num_objects>"
    exit()

# Peak height.
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

    # Generate locations
    x_vals = numpy.zeros(num_objects)
    y_vals = numpy.zeros(num_objects)
    z_vals = numpy.zeros(num_objects)
    h_vals = numpy.ones(num_objects) * intensity

    for j in range(num_objects):

        # Random locations.
        if 0:
            x_vals[j] = 0.05 * float(x_size) + 0.9 * float(x_size) * random.random()
            y_vals[j] = 0.05 * float(y_size) + 0.9 * float(y_size) * random.random()
            z_off = 0.8*(random.random() - 0.4)
            z_vals[j] = z_off * 1000.0

        # On a grid.
        if 1:
            #x_vals[j] = 10.0 + 10.0*(j%23)
            #y_vals[j] = 10.0 + 10.0*math.floor(float(j)/23.0)
            x_vals[j] = 20.0 + 20.0*(j%10)
            y_vals[j] = 20.0 + 20.0*math.floor(float(j)/10.0)            
            z_off = -0.5 + float(j)/float(num_objects - 1)
            #z_off = -0.4 + 0.8 * float(j)/float(num_objects - 1)
            #z_off = 0.0
            z_vals[j] = z_off * 1000.0

    # Generate objects.
    objects = PSF.PSF(x_vals, y_vals, z_vals, h_vals)
    
    # Draw the image.
    image = dg.drawGaussians([x_size, y_size], objects, background = 50, res = 5)
    #image = dg.drawGaussians([x_size, y_size], objects, background = 0, res = 5)
    #image[0:(x_size/2),:] += 50

    # Add poisson noise and baseline.
    image = numpy.random.poisson(image) + 100.0
        
    # Save the image.
    dax_data.addFrame(image)

    # Save the molecule locations.
    a_vals = PSF.PSFIntegral(z_vals, h_vals)
    
    ax = numpy.ones(num_objects)
    ww = 2.0 * 160.0 * numpy.ones(num_objects)
    if (PSF.psf_type == "astigmatic"):
        sx_vals = objects[:,3]
        sy_vals = objects[:,4]
        ax = sy_vals/sx_vals
        ww = 2.0*160.0*numpy.sqrt(sx_vals*sy_vals)
        
    mols = i3dtype.createDefaultI3Data(num_objects)
    i3dtype.posSet(mols, 'x', x_vals + 1.0)
    i3dtype.posSet(mols, 'y', y_vals + 1.0)
    i3dtype.posSet(mols, 'z', z_vals)
    i3dtype.setI3Field(mols, 'a', a_vals)
    i3dtype.setI3Field(mols, 'h', h_vals)
    i3dtype.setI3Field(mols, 'w', ww)
    i3dtype.setI3Field(mols, 'ax', ax)
    i3dtype.setI3Field(mols, 'fr', i+1)
    i3_data.addMolecules(mols)

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
