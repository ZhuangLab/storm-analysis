#!/usr/bin/env python
"""
sCMOS calibration tests.
"""
import numpy
import numpy.random
import tifffile

import storm_analysis

import storm_analysis.simulator.camera as camera
import storm_analysis.sCMOS.movie_to_calib_format as movieToCalFmt


def test_create_cal_1():

    cal_name = storm_analysis.getPathOutputTest("scmos.npy")

    x_size = 110
    y_size = 100
    gain = 2.0
    offset = 100.0
    read_noise = 1.0
    [cam_offset, cam_variance, cam_gain] = camera.createSCMOSCalibration(cal_name,
                                                                         x_size,
                                                                         y_size,
                                                                         gain,
                                                                         read_noise,
                                                                         hot_fraction = 0.0,
                                                                         hot_lambda = 10.0,
                                                                         offset = offset)
    
    # Check values.
    assert numpy.allclose(cam_offset, offset * numpy.ones((x_size, y_size)))
    assert numpy.allclose(cam_gain, gain * numpy.ones((x_size, y_size)))
    assert (abs(numpy.mean(cam_variance/(gain*gain)) - 1.0) < 0.1)

    
def test_create_cal_2():

    cal_name = storm_analysis.getPathOutputTest("scmos.npy")

    x_size = 110
    y_size = 100
    gain = 2.0
    offset = 100.0
    read_noise = 0.0
    [cam_offset, cam_variance, cam_gain] = camera.createSCMOSCalibration(cal_name,
                                                                         x_size,
                                                                         y_size,
                                                                         gain,
                                                                         read_noise,
                                                                         hot_fraction = 1.0,
                                                                         hot_lambda = 10.0,
                                                                         offset = offset)
    
    # Check values.
    assert numpy.allclose(cam_offset, offset * numpy.ones((x_size, y_size)))
    assert numpy.allclose(cam_gain, gain * numpy.ones((x_size, y_size)))
    assert (abs(numpy.mean(numpy.sqrt(cam_variance))/(gain * 10.0) - 1.0) < 0.1)

    
def test_mtcf_1():

    tif_name = storm_analysis.getPathOutputTest("mtcf.tif")

    rd_noise = numpy.ones((10,5))
    rd_noise[5:,:] = 2.0

    offset = 10.0*numpy.ones(rd_noise.shape)
    offset[3:,:] = 20.0

    # Create calibration movie.
    with tifffile.TiffWriter(tif_name) as tf:
        for i in range(1000):
            image = numpy.random.normal(scale = rd_noise, size = rd_noise.shape)
            image += offset
            tf.save(numpy.round(image).astype(numpy.uint16))

    [frame_mean, N, NN] = movieToCalFmt.movieToCalibration(tif_name)
    mean = N/float(frame_mean.size)
    variance = NN/float(frame_mean.size) - mean*mean
    rd_sqr = rd_noise * rd_noise
    
    assert(numpy.allclose(mean, offset, rtol = 0.1))
    assert(numpy.allclose(rd_sqr, variance, rtol = 0.5))
    

if (__name__ == "__main__"):
    test_create_cal_1()
    test_create_cal_2()
    test_mtcf_1()
    
    
