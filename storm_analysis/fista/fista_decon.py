#!/usr/bin/env python
"""
Deconvolve images in 3D using FISTA.

Hazen 1/16
"""
import numpy

import storm_analysis.fista.fista_3d as fista_3d
import storm_analysis.fista.fista_decon_utilities_c as fdUtil
import storm_analysis.fista.fista_fft_c as fistaFFTC

#
# FIXME: Ignore peaks outside of user-specified AOI and/or only
#        do decon on the specified sub-region.
#
class FISTADecon(object):

    def __init__(self, image_size, psf_object, number_zvals, timestep):
        """
        
        """
        self.background = numpy.zeros(image_size)
        self.new_image = None
        self.psf_object = psf_object

        # Calculate z values to use if 3D.
        if (self.psf_object.getType() == "3D"):
            self.z_min = self.psf_object.getZMin() + 1.0
            self.z_max = self.psf_object.getZMax() - 1.0
            z_step = (self.z_max - self.z_min)/float(number_zvals - 1.0)        
            self.zvals = []
            for i in range(number_zvals):
                self.zvals.append(self.z_min + float(i) * z_step)
        else:
            self.z_min = 0.0
            self.z_max = 1.0
            self.zvals = [0.0]

        # Create PSFs.
        size_x, size_y = image_size
        psfs = numpy.zeros((size_x, size_y, len(self.zvals)))
        for i in range(len(self.zvals)):
            psfs[:,:,i] = self.psf_object.getPSF(self.zvals[i],
                                                 shape = image_size,
                                                 normalize = True)

        if False:
            # Python solver (useful for debugging).
            print("Using Python solver.")
            self.fsolver = fista_3d.FISTA(psfs, timestep)
        else:
            # C solver (about 4x faster).
            print("Using C solver.")
            self.fsolver = fistaFFTC.FISTA(psfs, timestep)

    def cleanup(self):
        if isinstance(self.fsolver, fistaFFTC.FISTA):
            self.fsolver.cleanup()
        
    def decon(self, iterations, f_lambda, verbose = False):
        for i in range(iterations):
            if verbose and ((i%10) == 0):
                print(i, self.fsolver.l2Error())
            self.fsolver.iterate(f_lambda)

    def getPeaks(self, threshold, margin):
        """
        Extract peaks from the deconvolved image and create
        an array that can be used by a peak fitter.
        """
        fx = self.getXVector()

        fd_peaks = fdUtil.getPeaks(fx, threshold, margin)

        peaks = {"x" : fd_peaks[:,2],
                 "y" : fd_peaks[:,1]}

        if (fx.shape[2] > 1):

            # Initial z as the range 0.0 - 1.0.
            temp_z = fd_peaks[:,3]/(float(fx.shape[2])-1.0)

            # Convert z to nanometers.
            peaks["z"] = (self.z_max - self.z_min)*temp_z + self.z_min
            
        else:
            peaks["z"] = numpy.zeros(peaks["x"].size)
        
        return peaks
        
    def getXVector(self):
        return self.fsolver.getXVector()

    def getZRange(self):
        return [self.z_min, self.z_max]

    def newBackground(self, background):
        no_bg_image = self.new_image - background
        self.fsolver.newImage(no_bg_image)
        
    def newImage(self, image):
        self.new_image = image


#
# Deconvolution testing.
#
# FIXME: This is broken.
#
if (__name__ == "__main__"):

    import argparse
    import tifffile

    import storm_analysis.rolling_ball_bgr.rolling_ball as rollingBall
    import storm_analysis.sa_library.datareader as datareader
    import storm_analysis.sa_library.parameters as params
    import storm_analysis.sa_library.sa_h5py as saH5Py
    import storm_analysis.spliner.spline_to_psf as splineToPSF
    import storm_analysis.wavelet_bgr.wavelet_bgr as waveletBGR

    parser = argparse.ArgumentParser(description = 'FISTA deconvolution - Beck and Teboulle, SIAM J. Imaging Sciences, 2009')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie to deconvolve, can be .dax, .tiff or .spe format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")
    parser.add_argument('--output', dest='output', type=str, required=True,
                        help = "The name of the .tif file to save the results in.")

    args = parser.parse_args()

    # Load parameters
    parameters = params.ParametersSplinerFISTA().initFromFile(args.settings)

    # Open movie and load the first frame.
    movie_data = datareader.inferReader(args.movie)
    [x_size, y_size, z_size] = movie_data.filmSize()
    image = (movie_data.loadAFrame(0) - parameters.getAttr("camera_offset"))/parameters.getAttr("camera_gain")
    image = image.astype(numpy.float)

    # Load spline.
    psf_object = splineToPSF.loadSpline(parameters.getAttr("spline"))
    
    # Do FISTA deconvolution.
    fdecon = FISTADecon(image.shape,
                        psf_object,
                        parameters.getAttr("fista_number_z"),
                        parameters.getAttr("fista_timestep"))

    if False:
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
        
    fdecon.newImage(image)
    fdecon.newBackground(background)

    fdecon.decon(parameters.getAttr("fista_iterations"),
                 parameters.getAttr("fista_lambda"),
                 verbose = True)

    # Save results.
    fx = fdecon.getXVector()
    print(numpy.min(fx), numpy.max(fx))
    with tifffile.TiffWriter(args.output) as tf:
        tf.save(image.astype(numpy.float32))
        for i in range(fx.shape[2]):
            tf.save(fx[:,:,i].astype(numpy.float32))
    
    # Find peaks in the decon data.
    peaks = fdecon.getPeaks(parameters.getAttr("fista_threshold"), 5)

    saH5Py.saveLocalizations(args.output[:-4] + ".hdf5", peaks)


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
