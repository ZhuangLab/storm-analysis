#!/usr/bin/env python
"""
Rolling ball background estimation.

Hazen 02/16
"""

import numpy
import scipy
import scipy.ndimage

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.datawriter as datawriter
    
import storm_analysis.rolling_ball_bgr.rolling_ball_lib_c as rollingBallLibC
import storm_analysis.rolling_ball_bgr.rolling_ball_py as rollingBallPy


class RollingBall(rollingBallLibC.CRollingBall):
    """
    Rolling ball smoothing class.
    """
    pass


def rollingBallSub(movie_in, movie_out, radius, sigma, offset = 100):
        
    input_movie = datareader.inferReader(movie_in)
    output_dax = datawriter.inferWriter(movie_out)

    rb = RollingBall(radius, sigma)
        
    for i in range(input_movie.filmSize()[2]):

        if((i%10) == 0):
            print("Processing frame", i)

        image = input_movie.loadAFrame(i) - offset

        if False:
            image = image.astype(numpy.float)
            lowpass = scipy.ndimage.filters.gaussian_filter(image, sigma)
            sub = image - lowpass
            
        else:
            sub = rb.removeBG(image)
            
        output_dax.addFrame(sub + offset)

    output_dax.close()


if (__name__ == "__main__"):

    import argparse
    
    parser = argparse.ArgumentParser(description = 'Rolling ball background subtraction')

    parser.add_argument('--movie_in', dest='movie_in', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--movie_out', dest='movie_out', type=str, required=True,
                        help = "The name of the movie to save the results. This will always be .dax format.")
    parser.add_argument('--radius', dest='radius', type=float, required=True,
                        help = "The radius of the ball in pixels.")
    parser.add_argument('--sigma', dest='sigma', type=float, required=True,
                        help = "Sigma of the gaussian to use for image smoothing prior to the rolling ball step.")
    parser.add_argument('--baseline', dest='baseline', type=float, required=False, default=100,
                        help = "Camera baseline in ADU.")

    args = parser.parse_args()

    rollingBallSub(args.movie_in, args.movie_out, args.radius, args.sigma, args.baseline)
