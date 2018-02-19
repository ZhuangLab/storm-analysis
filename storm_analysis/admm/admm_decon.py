#!/usr/bin/python
#
# Uses ADMM to perform image deconvolution. It is assumed that
# the background has already been subtracted from the image.
#
# Hazen 02/15
#

import numpy

import admm_lasso_c

import sa_library.rebin as rebin

class ADMMDeconException(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

class ADMMDecon(object):

    ## __init__
    #
    # @param image_shape [xsize, ysize].
    # @param psf A 2D numpy array of size [xsize, ysize] or a floating point number that will be used as the sigma of a 2D symmetric gaussian.
    # @param rho The rho parameter.
    # @param a_lambda The lambda parameter.
    # @upscaling_factor An integer >= 1.
    #
    def __init__(self, image_shape, psf, rho, a_lambda, upscaling_factor):
        self.a_lambda = a_lambda
        self.rho = rho
        self.upscaling_factor = upscaling_factor

        if (image_shape[0] != image_shape[1]):
            raise ADMMDeconException("Images must be square!")

        self.xsize = image_shape[0] * upscaling_factor
        self.ysize = image_shape[1] * upscaling_factor

        # Check if this is a measured psf.
        if (type(psf) == type(numpy.array)):

            # Check measured psf size.
            if (psf.shape[0] != image_shape[0]) or (psf.shape[1] != image_shape[1]):
                raise ADMMDeconException("PSF size must match image size!")

            # Upscale to final size.
            self.psf = rebin.upSampleFFT(psf, upscaling_factor)

        # If not, then psf is assumed to be the sigma of a 2D-symmetric gaussian PSF.
        else:
            dx = numpy.arange(self.xsize) - (self.xsize)/2.0
            sigma = float(psf) * float(upscaling_factor)
            exp = numpy.exp(-1.0*dx*dx/(2.0*sigma*sigma))
            self.psf = numpy.outer(exp, exp)

        # Create lasso solver.
        self.lasso = admm_lasso_c.ADMMLasso(self.psf, self.rho)

    def cleanup(self):
        self.lasso.cleanup()

    def decon(self, iterations):
        self.lasso.iterate(self.a_lambda, iterations)
        return self.lasso.getXVector()

    def newImage(self, image):
        self.up_image = rebin.upSampleFFT(image, self.upscaling_factor) * self.upscaling_factor * self.upscaling_factor
        self.lasso.newImage(self.up_image)


# Testing
if (__name__ == "__main__"):

    import sys

    import sa_library.datareader as datareader
    import sa_library.daxwriter as daxwriter
    import wavelet_bgr.wavelet_bgr as wavelet_bgr

    dxr = datareader.inferReader(sys.argv[1])
    frame = dxr.loadAFrame(0) - 100

    wbgr = wavelet_bgr.WaveletBGR(wavelet_type = 'db6')
    dxw = daxwriter.DaxWriter("bg.dax", 0, 0)
    dxw.addFrame(frame)
    frame = wbgr.removeBG(frame, 5, 100, 5)
    dxw.addFrame(frame)
    dxw.close()

    a_lambda = 2.0 * 0.02 * numpy.max(frame)

    dxw = daxwriter.DaxWriter("test.dax", 0, 0)
    decon = ADMMDecon(frame.shape, 1.0, 0.05, a_lambda, 1)
    decon.newImage(frame)
    dxw.addFrame(decon.up_image)
    #dxw.addFrame(1000.0 * decon.psf)
    for i in range(10):
        result = decon.decon(10)
        dxw.addFrame(result)

    dxw.close()

    decon.cleanup()

#
# The MIT License
#
# Copyright (c) 2015 Zhuang Lab, Harvard University
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
