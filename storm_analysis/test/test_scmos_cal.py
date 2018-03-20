#!/usr/bin/env python
"""
sCMOS calibration tests.
"""
import numpy

import storm_analysis

import storm_analysis.simulator.camera as camera

import storm_analysis.test.verifications as veri


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
    

if (__name__ == "__main__"):
    test_create_cal_1()
    test_create_cal_2()
    
