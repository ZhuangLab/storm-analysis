#!/usr/bin/env python
#
# Rolling ball background estimation.
#
# Hazen 02/16
#

import numpy
import scipy
import scipy.ndimage

import storm_analysis.rolling_ball_bgr.rolling_ball_lib_c as rollingBallLibC
import storm_analysis.rolling_ball_bgr.rolling_ball_py as rollingBallPy

#
# Rolling ball smoothing class.
#
class RollingBall(rollingBallLibC.CRollingBall):
    pass

#class RollingBall(rollingBallPy.PyRollingBall):
#    pass

if (__name__ == "__main__"):

    import sys
    
    import storm_analysis.sa_library.datareader as datareader
    import storm_analysis.sa_library.daxwriter as daxwriter
        
    if (len(sys.argv) < 5):
        print("usage <movie (in)> <subtracted movie (out> <ball radius> <smoothing sigma> <baseline (optional, 100 default)>")
        exit()

    input_movie = datareader.inferReader(sys.argv[1])
    output_dax = daxwriter.DaxWriter(sys.argv[2], 0, 0)

    rb = RollingBall(float(sys.argv[3]), float(sys.argv[4]))

    offset = 100.0
    if (len(sys.argv) == 6):
        offset = float(sys.argv[5])
        
    for i in range(input_movie.filmSize()[2]):

        if((i%10) == 0):
            print("Processing frame", i)

        image = input_movie.loadAFrame(i) - offset

        if 0:
            image = image.astype(numpy.float)
            lowpass = scipy.ndimage.filters.gaussian_filter(image, float(sys.argv[3]))
            sub = image - lowpass
            
        else:
            sub = rb.removeBG(image)
            
        output_dax.addFrame(sub + offset)

    output_dax.close()
