#!/usr/bin/env python
#
# Deconvolve images in 3D using FISTA.
#
# Hazen 1/16
#

import numpy

import sa_library.ia_utilities_c as utilC
import sa_library.rebin as rebin
import simulator.drawgaussians as dg
import spliner.spline_to_psf as splineToPsf

import fista_3d
import fista_decon_utilities_c as fdUtil
    
class FISTADecon(object):

    #
    # Upsample is the multiplier to use for re-sizing the image,
    #    for example upsample = 2 means to enlarge by 2x.
    #
    def __init__(self, image_size, spline, zvals, upsample = 1, background_sigma = 10):
        self.psf_heights = []
        self.upsample = int(upsample)
        self.zvals = zvals
        
        im_size_x, im_size_y = image_size
        size_x = im_size_x * self.upsample
        size_y = im_size_y * self.upsample
        
        s_to_psf = splineToPsf.SplineToPSF(spline)
        spline_size_x = spline_size_y = s_to_psf.getSize()

        psfs = numpy.zeros((size_x, size_y, len(self.zvals) + 1))

        # Add background fitting PSF.
        objects = numpy.zeros((1,5))

        objects[0,:] = [float(size_x) * 0.5,
                        float(size_y) * 0.5,
                        1.0,
                        float(background_sigma * self.upsample),
                        float(background_sigma * self.upsample)]
        temp = dg.drawGaussians([size_x, size_y], objects, background = 0, res = 5)
        temp = temp/numpy.sum(temp)
        psfs[:,:,0] = temp
        
        # Add PSFs.
        start_x = im_size_x/2 - spline_size_x/4
        start_y = im_size_y/2 - spline_size_y/4
        end_x = start_x + spline_size_x
        end_y = start_y + spline_size_y        
        for i in range(len(self.zvals)):
            temp = numpy.zeros((im_size_x, im_size_y))
            temp[start_x:end_x,start_y:end_y] = s_to_psf.getPSF(zvals[i])

            if (self.upsample > 1):
                temp = rebin.upSampleFFT(temp, self.upsample)
                
            psfs[:,:,i+1] = temp/numpy.sum(temp)
            self.psf_heights.append(numpy.max(psfs[:,:,i+1]))

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

    # Extract peaks from the deconvolved image and create
    # an array that can be used by a peak fitter.
    def getPeaks(self, threshold, margin):
        fx = self.getX()

        # Get area, position, height.
        no_bg_fx = fx[:,:,1:]
        fd_peaks = fdUtil.getPeaks(no_bg_fx, threshold, margin)
        num_peaks = fd_peaks.shape[0]

        peaks = numpy.zeros((num_peaks, utilC.getNResultsPar()))

        peaks[:,utilC.getXWidthIndex()] = numpy.ones(num_peaks)
        peaks[:,utilC.getYWidthIndex()] = numpy.ones(num_peaks)
        
        peaks[:,utilC.getXCenterIndex()] = fd_peaks[:,2]
        peaks[:,utilC.getYCenterIndex()] = fd_peaks[:,1]

        # Calculate height.
        h_index = utilC.getHeightIndex()
        for i in range(num_peaks):
            peaks[i,h_index] = fd_peaks[i,0] *self.psf_heights[int(round(fd_peaks[i,3]))]

        # Calculate z.
        peaks[:,utilC.getZCenterIndex()] = 1.0e-3 * ((fd_peaks[:,3]/(float(no_bg_fx.shape[2])-1.0)) * (self.zvals[-1] - self.zvals[0]) + self.zvals[0])

        # Background term calculation.
        # ..
        
        return peaks
        
    def getX(self):
        return self.fsolver.getX()

    def newImage(self, image):
        f_lambda = 20.0
        f_timestep = 1.0e-1
        if (self.upsample > 1):
            image = rebin.upSampleFFT(image, self.upsample)
        self.fsolver.newImage(image, f_lambda, f_timestep)
        

#
# Deconvolution testing.
#

if (__name__ == "__main__"):

    import sys
    
    import sa_library.datareader as datareader
    import sa_library.daxwriter as daxwriter
    #import sa_library.i3dtype as i3dtype    
    import sa_library.writeinsight3 as writeinsight3
    
    import fista_decon_utilities_c as fdUtil
    
    if (len(sys.argv) != 5):
        print "usage: <movie, input> <spline, input> <epsilon, input> <decon, output>"
        exit()

    # The Z planes to use.
    z_values = [-500.0, -250.0, 0.0, 250.0, 500.0]

    movie_data = datareader.inferReader(sys.argv[1])
    [x_size, y_size, z_size] = movie_data.filmSize()
    epsilon = float(sys.argv[3])

    image = movie_data.loadAFrame(0) - 100

    # Do FISTA deconvolution.
    fdecon = FISTADecon(image.shape, sys.argv[2], z_values, upsample = 1)
    fdecon.newImage(image)
    fdecon.decon()

    # Save results.
    fx = fdecon.getX()
    decon_data = daxwriter.DaxWriter(sys.argv[4], fx.shape[0], fx.shape[1])
    for i in range(fx.shape[2]):
        decon_data.addFrame(fx[:,:,i])
    decon_data.close()
    
    # Find peaks in the decon data.
    peaks = fdecon.getPeaks(1.0, 0)

    i3_writer = writeinsight3.I3Writer(sys.argv[4][:-4] + "_flist.bin")    
    #mols = i3dtype.createDefaultI3Data(px.size)
    #i3dtype.posSet(mols, 'x', px + 1.0)
    #i3dtype.posSet(mols, 'y', py + 1.0)
    #i3dtype.posSet(mols, 'z', pz)
    #i3dtype.setI3Field(mols, 'a', pi)
    #i3dtype.setI3Field(mols, 'h', ph)
    #i3dtype.setI3Field(mols, 'bg', pb)
    #i3_writer.addMolecules(mols)
    i3_writer.addMultiFitMolecules(peaks, x_size, y_size, 1, 160.0)
    i3_writer.close()


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
