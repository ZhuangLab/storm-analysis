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
import pickle
import tifffile

import storm_analysis.multi_plane.mp_utilities as mpUtil

import storm_analysis.sa_library.affine_transform_c as affineTransformC
import storm_analysis.sa_library.analysis_io as analysisIO
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
        bg_frame = frame.copy()
        bg_frame[mask] = bg_estimate[mask]

        # Update foreground and background estimates.
        bg_estimate = bg_filter.convolve(bg_frame)
        fg_estimate = fg_filter.convolve(frame - bg_estimate)        

    return bg_estimate


def loadImage(movie, frame, offset, gain, transform = None):
    """
    Load a single frame, apply sCMOS correction and remove
    values that are less than 1.0.
    """
    image = movie.loadAFrame(frame)
    image = (image - offset)*gain

    # Apply transform if requested.
    if transform is not None:
        image = transform.transform(image)

    # Remove values that are less than 1.0 as these will
    # trip up the square root function.
    mask = (image < 1.0)
    if (numpy.sum(mask) > 0):
        image[mask] = 1.0

    return image
    

def findOffsets(base_name, params_file, background_scale = 4.0, foreground_scale = 1.0, im_slice = None):
    """
    The 'main' function of this module.

    base_name - The basename for the group of movies.
    params_file - An analysis XML file containing the details for this experiment.
    background_scale - Features in the background change on this scale (in pixels)
                       or more slowly.
    foreground_scale - Features that change on this scale are likely foreground.
    im_slice - A slice object created for example with numpy.s_ to limit the analysis
               to a smaller AOI.

    Notes: 
      1. This only checks a limited range of offsets between the two channels.
      2. This assumes that the movies are longer than just a few frames.
    """
    n_tests = 10
    search_range = 5
    
    # Load parameters.
    parameters = params.ParametersMultiplane().initFromFile(params_file)

    # Load the movies from each camera.
    n_channels = 0
    movies = []
    for ext in mpUtil.getExtAttrs(parameters):
        movie_name = base_name + parameters.getAttr(ext)
        movies.append(datareader.inferReader(movie_name))
        n_channels += 1

    print("Found", n_channels, "movies.")

    # Load sCMOS calibration data.
    offsets = []
    gains = []
    for calib_name in mpUtil.getCalibrationAttrs(parameters):
        [offset, variance, gain, rqe] = analysisIO.loadCMOSCalibration(parameters.getAttr(calib_name))
        offsets.append(offset)
        gains.append(1.0/gain)

    assert(len(offsets) == n_channels)
    
    # Load the plane to plane mapping data & create affine transform objects.
    mappings = {}
    with open(parameters.getAttr("mapping"), 'rb') as fp:
        mappings = pickle.load(fp)

    atrans = []
    for i in range(n_channels-1):
        xt = mappings["0_" + str(i+1) + "_x"]
        yt = mappings["0_" + str(i+1) + "_y"]
        atrans.append(affineTransformC.AffineTransform(xt = xt, yt = yt))

    # Create background and foreground variance filters.
    #
    # FIXME: Is this right for movies that are not square?
    #
    [y_size, x_size] = movies[0].filmSize()[:2]

    if im_slice is not None:
        y_size = im_slice[0].stop - im_slice[0].start
        x_size = im_slice[1].stop - im_slice[1].start
        
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
    if False:
        frame = loadImage(movies[0], 0, offsets[0], gain[0])
        frame_bg = estimateBackground(frame, bg_filter, fg_filter, var_filter)
        with tifffile.TiffWriter("bg_estimate.tif") as tif:
            tif.save(frame.astype(numpy.float32))
            tif.save(frame_bg.astype(numpy.float32))
            tif.save((frame - frame_bg).astype(numpy.float32))

    votes = numpy.zeros((n_channels - 1, 2*search_range+1))
    for i in range(n_tests):
        print("Test", i)
        
        # Load reference frame.
        ref_frame = loadImage(movies[0], search_range + i, offsets[0], gain[0])
        if im_slice is not None:
            ref_frame = ref_frame[im_slice]
        ref_frame_bg = estimateBackground(ref_frame, bg_filter, fg_filter, var_filter)
        ref_frame -= ref_frame_bg

        # Load test frames and measure correlation.
        for j in range(n_channels - 1):
            best_corr = 0.0
            best_offset = 0
            for k in range(-search_range, search_range + 1):
                test_frame = loadImage(movies[j+1], search_range + i + k, offsets[j+1], gain[j+1], transform = atrans[j])
                if im_slice is not None:
                    test_frame = test_frame[im_slice]
                test_frame_bg = estimateBackground(test_frame, bg_filter, fg_filter, var_filter)
                test_frame -= test_frame_bg
                test_frame_corr = numpy.sum(ref_frame*test_frame)/numpy.sum(test_frame)
                if (test_frame_corr > best_corr):
                    best_corr = test_frame_corr
                    best_offset = k + search_range

            votes[j, best_offset] += 1

    # Print results.
    print("Offset votes:")
    print(votes)

    frame_offsets = [0]
    frame_offsets += list(numpy.argmax(votes, axis = 1) - search_range)
    print("Best offsets:")
    for i in range(n_channels):
        print(str(i) + ": " + str(frame_offsets[i]))

    # Create stacks with optimal offsets.
    print("Saving image stacks.")
    for i in range(n_channels):
        with tifffile.TiffWriter(base_name + "_offsets_ch" + str(i) + ".tif") as tif:
            for j in range(5):
                if (i == 0):
                    frame = loadImage(movies[i],
                                      search_range + frame_offsets[i] + j,
                                      offsets[i],
                                      gain[i])
                else:
                    frame = loadImage(movies[i],
                                      search_range + frame_offsets[i] + j,
                                      offsets[i],
                                      gain[i],
                                      transform = atrans[i-1])
                if im_slice is not None:
                    frame = frame[im_slice]
                    
                frame_bg = estimateBackground(frame, bg_filter, fg_filter, var_filter)
                frame -= frame_bg
                tif.save(frame.astype(numpy.float32))

    return frame_offsets
    

if (__name__ == "__main__"):

    import argparse    

    parser = argparse.ArgumentParser(description = 'Estimate offsets between movies.')

    parser.add_argument('--basename', dest='basename', type=str, required=True,
                        help = "The base name of the movie to analyze.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    findOffsets(args.basename, args.settings)
    
