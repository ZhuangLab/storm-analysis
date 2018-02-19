#!/usr/bin/env python
"""
Uses ADMM to perform image deconvolution. 

Hazen 02/18
"""
import numpy

import storm_analysis.sa_library.cs_decon as csDecon

import storm_analysis.admm.admm_lasso_c as admmLassoC


class ADMMDecon(csDecon.CSDecon):
    
    def __init__(self, image_size, psf_object, number_zvals, rho):
        super(ADMMDecon, self).__init__(image_size, psf_object, number_zvals)

        assert(number_zvals == 1), "Only 2D deconvolution is currently supported."

        psfs = self.createPSFs()
        
        self.cs_solver = admmLassoC.ADMMLasso(psfs, rho)


# Deconvolution testing. For now we are just copying the FISTA equivalent. We're
# also using the FISTA parameters instead of creating a new ADMM equivalent set.
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

    parser = argparse.ArgumentParser(description = 'ADMM deconvolution - Boyd et al, Foundations and Trends in Machine Learning, 2011.')

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
    adecon = ADMMDecon(image.shape,
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
        
    adecon.newImage(image)
    adecon.newBackground(background)

    adecon.decon(parameters.getAttr("fista_iterations"),
                 parameters.getAttr("fista_lambda"),
                 verbose = True)

    # Save results.
    fx = adecon.getXVector()
    print("x vector min, max", numpy.min(fx), numpy.max(fx))
    with tifffile.TiffWriter(args.output) as tf:
        tf.save(image.astype(numpy.float32))
        for i in range(fx.shape[2]):
            tf.save(fx[:,:,i].astype(numpy.float32))
            
    # Find peaks in the decon data.
    peaks = adecon.getPeaks(parameters.getAttr("fista_threshold"), 5)

    saH5Py.saveLocalizations(args.output[:-4] + ".hdf5", peaks)

    # Clean up
    adecon.cleanup()

#
# The MIT License
#
# Copyright (c) 2018 Babcock Lab, Harvard University
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
