#!/usr/bin/env python
"""
In multiple camera setups getting them all to start at
the same time can be tricky. This attempts to automatically
determine the frame number offset between the movies
from different cameras. It uses a correlation based
approach, so there should be at least a few features in 
common in the images from the different cameras.

Note: This assumes Poisson statistics when estimating
      the background.

Hazen 08/17
"""
import numpy
import tifffile

import storm_analysis.multi_plane.mp_utilities_c as mpUtilC

import storm_analysis.sa_library.affine_transform_c as affineTransformC
import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.matched_filter_c as matchedFilterC
import storm_analysis.sa_library.parameters as params

import storm_analysis.simulator.draw_gaussians_c as dg

import storm_analysis.spliner.spline_to_psf as splineToPSF


def estimateBackground(frame, bg_filter, fg_filter, var_filter, reps = 5, threshold = 2.0):
    """
    frame - numpy array containing the image to perform background estimation on.
    bg_filter - MatchedFilter object for estimating the background.
    fg_filter - MatchedFilter object for estimating the foreground.
    var_filter - MatchedFilter object for calculating the variance given the
                 smoothing of fg_filter.
    reps - Number of iterations of estimation improvement.
    threshold - Foreground versus background significance in units of sigma.
    """

    # Smooth image as an initial estimate of the background.
    bg_estimate = bg_filter.convolve(frame)

    # Smooth background subtracted image as initial
    # estimation of the foreground.
    fg_estimate = fg_filter.convolve(frame - bg_estimate)

    # Iterative update of the background estimate.
    for i in range(reps):

        # Calculate variance of the background estimate.
        bg_variance = var_filter.convolve(bg_estimate)

        #
        # Identify regions of the image that are likely foreground
        # based on the fact that their intensity is 'significant'
        # i.e. threshold sigma above the background estimate.
        #
        fg_bg_ratio = fg_estimate/numpy.sqrt(bg_variance)
        mask = (fg_bg_ratio > threshold)

        # Replace these regions with the current background estimate.
        bg_frame = frame[~mask] + bg_estimate[mask]

        # Update foreground and background estimates.
        bg_estimate = bg_filter.convolve(bg_frame)
        fg_estimate = fg_filter.convolve(frame - bg_estimate)        

    return bg_estimate

def loadImage(movie, frame, offset, gain):
    """
    Load a single frame, apply sCMOS correction and remove
    values that are less than 1.0.
    """
    image = movie.getFrame(frame)
    image = (image - offset)*gain
    mask = (image < 1.0)
    if (numpy.sum(mask) > 0):
        image[mask] = 1.0
    return image
    

def findOffsets(base_name, params_file, background_scale = 8.0, foreground_scale = 1.0):

    # Load parameters.
    parameters = params.ParametersMultiplane().initFromFile(params_file)

    # Load the movies from each camera.
    n_channels = 0
    movies = []
    for ext in mpUtilC.getExtAttrs(parameters):
        movie_name = base_name + parameters.getAttr(ext)
        movies.append(datareader.inferReader(movie_name))
        n_channels += 1

    print("Found", n_channels, "movies.")

    # Load sCMOS calibration data.
    offsets = []
    gains = []
    for calib_name in mpUtilC.getCalibrationAttrs(parameters):
        [offset, variance, gain] = numpy.load(calib_name)
        offsets.append(offset)
        gains.append(1.0/gain)

    assert(len(offsets) == n_channels)
    
    # Load the plane to plane mapping data & create affine transform objects.
    mappings = {}
    with open(parameters.getAttr("mapping"), 'rb') as fp:
        mappings = pickle.load(fp)

    # Subtract 1, because we added 1 to the x,y coordinates when we saved them.
    atrans = []
    for i in range(n_channels-1):
        xt = mpUtilC.marginCorrect(mappings["0_" + str(i+1) + "_x"], -1)
        yt = mpUtilC.marginCorrect(mappings["0_" + str(i+1) + "_y"], -1)
        atrans.append(affineTransformC.AffineTransform(xt = xt, yt = yt))

    # Create background and foreground variance filters.
    #
    # FIXME: Is this right for movies that are not square?
    #
    [y_size, x_size] = movies[0].getFilmSize()[:2]

    psf = dg.drawGaussiansXY((x_size, y_size),
                             numpy.array([0.5*x_size]),
                             numpy.array([0.5*y_size]),
                             sigma = background_scale)
    psf = psf/numpy.sum(psf)
    bg_filter = matchedFilterC.MatchedFilter(psf)

    psf = dg.drawGaussiansXY((x_size, y_size),
                             numpy.array([0.5*x_size]),
                             numpy.array([0.5*y_size]),
                             sigma = foreground_scale)
    psf = psf/numpy.sum(psf)
    fg_filter = matchedFilterC.MatchedFilter(psf)
    var_filter = matchedFilterC.MatchedFilter(psf*psf)

    # Check background estimation.
    if True:
        frame = loadImage(movies[0], 0, offsets[0], gain[0])
        frame_bg = estimateBackground(frame, bg_filter, fg_filter, var_filter)
        with tifffile.TiffWriter("bg_estimate.tif") as tif:
            tif.save(frame.astype(numpy.float32))
            tif.save(frame_bg.astype(numpy.float32))
            

if (__name__ == "__main__"):

    import argparse    

    parser = argparse.ArgumentParser(description = 'Estimate offsets between movies.')

    parser.add_argument('--basename', dest='basename', type=str, required=True,
                        help = "The base name of the movie to analyze.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    findOffsets(args.basename, args.settings)
    
