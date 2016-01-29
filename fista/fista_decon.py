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

#
# FIXME: Handling of non-square images.
#
# FIXME: Ignore peaks outside of user-specified AOI and/or only
#        do decon on the specified sub-region.
#
class FISTADecon(object):

    #
    # Upsample is the multiplier to use for re-sizing the image,
    #    for example upsample = 2 means to enlarge by 2x.
    #
    def __init__(self, image_size, spline_file, number_zvals, upsample = 1):
        self.background = numpy.zeros(image_size)
        self.psf_heights = []
        self.upsample = int(upsample)
        
        im_size_x, im_size_y = image_size
        size_x = im_size_x * self.upsample
        size_y = im_size_y * self.upsample
        
        s_to_psf = splineToPsf.SplineToPSF(spline_file)
        spline_size_x = spline_size_y = s_to_psf.getSize()

        # Calculate z values to use.
        z_min = s_to_psf.getZMin()
        z_max = s_to_psf.getZMax()
        z_step = (z_max - z_min)/float(number_zvals - 1.0)        
        self.zvals = []
        for i in range(number_zvals):
            self.zvals.append(z_min + float(i) * z_step)

        psfs = numpy.zeros((size_x, size_y, len(self.zvals)))

        # Add PSFs.
        start_x = im_size_x/2 - spline_size_x/4
        start_y = im_size_y/2 - spline_size_y/4
        end_x = start_x + spline_size_x/2
        end_y = start_y + spline_size_y/2
        for i in range(len(self.zvals)):
            temp = numpy.zeros((im_size_x, im_size_y))
            temp[start_x:end_x,start_y:end_y] = s_to_psf.getPSF(self.zvals[i])

            if (self.upsample > 1):
                temp = rebin.upSampleFFT(temp, self.upsample)
                
            psfs[:,:,i] = temp/numpy.sum(temp)
            self.psf_heights.append(numpy.max(psfs[:,:,i]))

        # Check PSFs.
        if 1:
            import sa_library.daxwriter as daxwriter

            psf_data = daxwriter.DaxWriter("fista_decon_psf.dax", psfs.shape[0], psfs.shape[1])
            for i in range(psfs.shape[2]):
                temp = psfs[:,:,i]
                psf_data.addFrame(1000.0 * temp/numpy.max(temp))
            psf_data.close()
            
        self.fsolver = fista_3d.FISTA(psfs)
            
    def decon(self, iterations, verbose = True):
        for i in range(iterations):
            if verbose and ((i%10) == 0):
                print i, self.fsolver.l2error()
            self.fsolver.iterate()

    # Extract peaks from the deconvolved image and create
    # an array that can be used by a peak fitter.
    #
    # FIXME: Margin should be 1/2 the spline size?
    #
    # FIXME: Need to compensate for up-sampling parameter in x,y.
    #
    def getPeaks(self, threshold, margin):
        fx = self.getX()

        # Get area, position, height.
        fd_peaks = fdUtil.getPeaks(fx, threshold, margin)
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
        peaks[:,utilC.getZCenterIndex()] = 1.0e-3 * ((fd_peaks[:,3]/(float(fx.shape[2])-1.0)) * (self.zvals[-1] - self.zvals[0]) + self.zvals[0])

        # Background term calculation.
        bg_index = utilC.getBackgroundIndex()
        for i in range(num_peaks):
            ix = int(round(fd_peaks[i,1]))
            iy = int(round(fd_peaks[i,2]))
            peaks[i,bg_index] = self.background[ix, iy]
            
        return peaks
        
    def getX(self):
        return self.fsolver.getX()

    def newImage(self, image, background, f_lambda, timestep):
        self.background = background

        no_bg_image = image - self.background
        if (self.upsample > 1):
            no_bg_image = rebin.upSampleFFT(no_bg_image, self.upsample)
        self.fsolver.newImage(no_bg_image, f_lambda, timestep)


#
# Deconvolution testing.
#
# FIXME: Use a parameters file to make user testing easier.
#
if (__name__ == "__main__"):

    import sys

    import sa_library.datareader as datareader
    import sa_library.daxwriter as daxwriter
    import sa_library.parameters as params
    import sa_library.writeinsight3 as writeinsight3
    import wavelet_bgr.wavelet_bgr as waveletBGR
    
    import fista_decon_utilities_c as fdUtil
    
    if (len(sys.argv) != 4):
        print "usage: <movie, input> <parameters_file, input> <decon, output>"
        exit()

    # Load parameters
    parameters = params.Parameters(sys.argv[2])

    # Open movie and load the first frame.
    movie_data = datareader.inferReader(sys.argv[1])
    [x_size, y_size, z_size] = movie_data.filmSize()
    image = movie_data.loadAFrame(0) - parameters.baseline
    image = image.astype(numpy.float)

    # Do FISTA deconvolution.
    fdecon = FISTADecon(image.shape,
                        parameters.spline,
                        parameters.fista_number_z,
                        upsample = parameters.fista_upsample)
    wbgr = waveletBGR.WaveletBGR()
    background = wbgr.estimateBG(image,
                                 parameters.wbgr_iterations,
                                 parameters.wbgr_threshold,
                                 parameters.wbgr_wavelet_level)
    fdecon.newImage(image,
                    background,
                    parameters.fista_lambda,
                    parameters.fista_timestep)
    fdecon.decon(parameters.fista_iterations)

    # Save results.
    fx = fdecon.getX()
    decon_data = daxwriter.DaxWriter(sys.argv[3], fx.shape[0], fx.shape[1])
    for i in range(fx.shape[2]):
        decon_data.addFrame(fx[:,:,i])
    decon_data.close()
    
    # Find peaks in the decon data.
    peaks = fdecon.getPeaks(parameters.fista_threshold, 0)

    i3_writer = writeinsight3.I3Writer(sys.argv[3][:-4] + "_flist.bin")    
    i3_writer.addMultiFitMolecules(peaks, x_size, y_size, 1, parameters.pixel_size)
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
