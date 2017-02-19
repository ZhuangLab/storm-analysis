#!/usr/bin/env python
"""
Deconvolve images in 3D using FISTA.

Hazen 1/16
"""

import numpy

import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.simulator.draw_gaussians_c as dg
import storm_analysis.spliner.spline_to_psf as splineToPSF

import storm_analysis.fista.fista_3d as fista_3d
import storm_analysis.fista.fista_decon_utilities_c as fdUtil
import storm_analysis.fista.fista_fft_c as fistaFFTC

#
# FIXME: Ignore peaks outside of user-specified AOI and/or only
#        do decon on the specified sub-region.
#
class FISTADecon(object):

    def __init__(self, image_size, spline_file, number_zvals, timestep, upsample = 1, check_psf = True):
        """
        Upsample is the multiplier to use for re-sizing the image,
        for example upsample = 2 means to enlarge by 2x.
        """
        
        self.background = numpy.zeros(image_size)
        self.psf_heights = []
        self.upsample = int(upsample)
        
        im_size_x, im_size_y = image_size
        size_x = im_size_x * self.upsample
        size_y = im_size_y * self.upsample

        # Load spline.
        s_to_psf = splineToPSF.loadSpline(spline_file)

        # Get size in X and Y.
        self.spline_size_x = self.spline_size_y = s_to_psf.getSize()

        # Calculate z values to use if 3D.
        if (s_to_psf.getType() == "3D"):
            self.z_min = s_to_psf.getZMin()
            self.z_max = s_to_psf.getZMax()
            z_step = (self.z_max - self.z_min)/float(number_zvals - 1.0)        
            self.zvals = []
            for i in range(number_zvals):
                self.zvals.append(self.z_min + float(i) * z_step)
        else:
            self.z_min = 0.0
            self.z_max = 1.0
            self.zvals = [0.0]

        psfs = numpy.zeros((size_x, size_y, len(self.zvals)))

        # Add PSFs.
        for i in range(len(self.zvals)):
            psfs[:,:,i] = s_to_psf.getPSF(self.zvals[i],
                                          shape = (im_size_x, im_size_y),
                                          up_sample = upsample,
                                          normalize = True)
            self.psf_heights.append(numpy.max(psfs[:,:,i]))
            #print "fista_decon", i, numpy.max(psfs[:,:,i])

        # Check PSFs.
        if check_psf:
            import os            
            import storm_analysis.sa_library.daxwriter as daxwriter

            psf_data = daxwriter.DaxWriter(os.path.join(os.path.dirname(spline_file), "fista_decon_psf.dax"),
                                           psfs.shape[0],
                                           psfs.shape[1])
            for i in range(psfs.shape[2]):
                temp = psfs[:,:,i]
                psf_data.addFrame(1000.0 * temp/numpy.max(temp))
            psf_data.close()

        if 0:
            # Python solver (useful for debugging).
            print("Using Python solver.")
            self.fsolver = fista_3d.FISTA(psfs, timestep)
        else:
            # C solver (about 4x faster).
            print("Using C solver.")
            self.fsolver = fistaFFTC.FISTA(psfs, timestep)

    def decon(self, iterations, f_lambda, verbose = False):
        for i in range(iterations):
            if verbose and ((i%10) == 0):
                print(i, self.fsolver.l2Error())
            self.fsolver.iterate(f_lambda)

    def getPeaks(self, threshold, margin):
        """
        Extract peaks from the deconvolved image and create
        an array that can be used by a peak fitter.
        
        FIXME: Need to compensate for up-sampling parameter in x,y.
        """
        
        fx = self.getXVector()

        # Get area, position, height.
        fd_peaks = fdUtil.getPeaks(fx, threshold, margin)
        num_peaks = fd_peaks.shape[0]

        peaks = numpy.zeros((num_peaks, utilC.getNPeakPar()))

        peaks[:,utilC.getXWidthIndex()] = numpy.ones(num_peaks)
        peaks[:,utilC.getYWidthIndex()] = numpy.ones(num_peaks)
        
        peaks[:,utilC.getXCenterIndex()] = fd_peaks[:,2]
        peaks[:,utilC.getYCenterIndex()] = fd_peaks[:,1]

        # Calculate height.
        #
        # FIXME: Typically the starting value for the peak height will be
        #        under-estimated unless a large enough number of FISTA
        #        iterations is performed to completely de-convolve the image.
        #
        h_index = utilC.getHeightIndex()
        #peaks[:,h_index] = fd_peaks[:,0]
        for i in range(num_peaks):
            peaks[i,h_index] = fd_peaks[i,0] * self.psf_heights[int(round(fd_peaks[i,3]))]

        # Calculate z (0.0 - 1.0).
        if (fx.shape[2] > 1):
            peaks[:,utilC.getZCenterIndex()] = fd_peaks[:,3]/(float(fx.shape[2])-1.0)

        # Background term calculation.
        bg_index = utilC.getBackgroundIndex()
        for i in range(num_peaks):
            ix = int(round(fd_peaks[i,1]))
            iy = int(round(fd_peaks[i,2]))
            peaks[i,bg_index] = self.background[ix, iy]
            
        return peaks
        
    def getXVector(self):
        return self.fsolver.getXVector()

    def getZRange(self):
        return [self.z_min, self.z_max]

    def newImage(self, image, background):
        self.background = background

        no_bg_image = image - self.background
        if (self.upsample > 1):
            no_bg_image = rebin.upSampleFFT(no_bg_image, self.upsample)
        self.fsolver.newImage(no_bg_image)


#
# Deconvolution testing.
#
if (__name__ == "__main__"):

    import sys

    import storm_analysis.rolling_ball_bgr.rolling_ball as rollingBall
    import storm_analysis.sa_library.datareader as datareader
    import storm_analysis.sa_library.daxwriter as daxwriter
    import storm_analysis.sa_library.parameters as params
    import storm_analysis.sa_library.writeinsight3 as writeinsight3
    import storm_analysis.wavelet_bgr.wavelet_bgr as waveletBGR
    
    import fista_decon_utilities_c as fdUtil
    
    if (len(sys.argv) != 4):
        print("usage: <movie, input> <parameters_file, input> <decon, output>")
        exit()

    # Load parameters
    parameters = params.ParametersSplinerFISTA().initFromFile(sys.argv[2])

    # Open movie and load the first frame.
    movie_data = datareader.inferReader(sys.argv[1])
    [x_size, y_size, z_size] = movie_data.filmSize()
    image = movie_data.loadAFrame(0) - parameters.getAttr("baseline")
    image = image.astype(numpy.float)

    # Do FISTA deconvolution.
    fdecon = FISTADecon(image.shape,
                        parameters.getAttr("spline"),
                        parameters.getAttr("fista_number_z"),
                        parameters.getAttr("fista_timestep"),
                        upsample = parameters.getAttr("fista_upsample"))

    if 0:
        # Wavelet background removal.
        wbgr = waveletBGR.WaveletBGR()
        background = wbgr.estimateBG(image,
                                     parameters.getAttr("wbgr_iterations"),
                                     parameters.getAttr("wbgr_threshold"),
                                     parameters.getAttr("wbgr_wavelet_level"))
    else:
        # Rolling ball background removal.
        rb = rollingBall.RollingBall(parameters.getAttr("rb_radius"),
                                     parameters.getAttr("rb_sigma"))
        background = rb.estimateBG(image)
        
    fdecon.newImage(image, background)

    fdecon.decon(parameters.getAttr("fista_iterations"),
                 parameters.getAttr("fista_lambda"),
                 verbose = True)

    # Save results.
    fx = fdecon.getXVector()
    print(numpy.min(fx), numpy.max(fx))
    decon_data = daxwriter.DaxWriter(sys.argv[3], fx.shape[0], fx.shape[1])
    for i in range(fx.shape[2]):
        decon_data.addFrame(fx[:,:,i])
    decon_data.close()
    
    # Find peaks in the decon data.
    peaks = fdecon.getPeaks(parameters.getAttr("threshold"), 5)

    zci = utilC.getZCenterIndex()
    z_min, z_max = fdecon.getZRange()
    peaks[:,zci] = 1.0e-3 * ((z_max - z_min)*peaks[:,zci] + z_min)
    
    i3_writer = writeinsight3.I3Writer(sys.argv[3][:-4] + "_flist.bin")    
    i3_writer.addMultiFitMolecules(peaks, x_size, y_size, 1, parameters.getAttr("pixel_size"))
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
