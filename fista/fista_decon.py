#!/usr/bin/env python
#
# Deconvolve images in 3D using FISTA.
#
# Hazen 1/16
#

import numpy

import sa_library.rebin as rebin
import spliner.spline_to_psf as spline_to_psf

import fista_3d

class FISTADecon(object):

    #
    # Upsample is the multiplier to use for re-sizing the image,
    #    for example upsample = 2 means to enlarge by 2x.
    #
    def __init__(self, image_size, spline, zvals, upsample = 1):
        self.upsample = int(upsample)
        
        im_size_x, im_size_y = image_size
        
        s_to_psf = spline_to_psf.SplineToPSF(spline)
        spline_size_x = spline_size_y = s_to_psf.getSize()

        start_x = im_size_x/2 - spline_size_x/4
        start_y = im_size_y/2 - spline_size_y/4
        end_x = start_x + spline_size_x
        end_y = start_y + spline_size_y

        psfs = numpy.zeros((im_size_x * self.upsample, im_size_y * self.upsample, len(zvals)))
        for i in range(len(zvals)):
            temp = numpy.zeros((im_size_x, im_size_y))
            temp[start_x:end_x,start_y:end_y] = s_to_psf.getPSF(zvals[i])

            if (self.upsample > 1):
                temp = rebin.upSampleFFT(temp, self.upsample)
                
            psfs[:,:,i] = temp/numpy.sum(temp)

        # Check PSFs.
        if 1:
            import sa_library.daxwriter as daxwriter

            psf_data = daxwriter.DaxWriter("fista_decon_psf.dax", psfs.shape[0], psfs.shape[1])
            for i in range(psfs.shape[2]):
                temp = psfs[:,:,i]
                psf_data.addFrame(1000.0 * temp/numpy.max(temp))
            psf_data.close()
            
        self.fsolver = fista_3d.FISTA(psfs)
            
    def decon(self, iterations = 100, verbose = True):
        for i in range(iterations):
            if verbose and ((i%10) == 0):
                print i, self.fsolver.l2error()
            self.fsolver.iterate()

    def getX(self):
        return self.fsolver.getX()

    def newImage(self, image):
        f_lambda = 20.0
        f_timestep = 1.0e-1
        if (self.upsample > 1):
            image = rebin.upSampleFFT(image, self.upsample)
        self.fsolver.newImage(image, f_lambda, f_timestep)
        

# Deconvolution testing.
if (__name__ == "__main__"):

    import pickle
    import sys
    
    import sa_library.datareader as datareader
    import sa_library.daxwriter as daxwriter

    if (len(sys.argv) != 5):
        print "usage: <movie, input> <spline, input> <epsilon, input> <decon, output>"
        exit()

    # The Z planes to use.
    z_values = [-500.0, -250.0, 0.0, 250.0, 500.0]

    movie_data = datareader.inferReader(sys.argv[1])
    epsilon = float(sys.argv[3])

    image = movie_data.loadAFrame(0) - 100

    # Do FISTA deconvolution.
    fdecon = FISTADecon(image.shape, sys.argv[2], z_values, upsample = 1)
    fdecon.newImage(image)
    fdecon.decon()

    # Save results.
    fx = fdecon.getX()
    fx = fx/numpy.max(fx)
    decon_data = daxwriter.DaxWriter(sys.argv[4], fx.shape[0], fx.shape[1])
    for i in range(fx.shape[2]):
        decon_data.addFrame(1000.0 * fx[:,:,i])
    decon_data.close()

    
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
