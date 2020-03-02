#!/usr/bin/env python
"""
sCMOS calibration tests.
"""
import numpy
import numpy.random
import tifffile

import storm_analysis

import storm_analysis.simulator.camera as camera
import storm_analysis.sCMOS.camera_calibration as camCal
import storm_analysis.sCMOS.movie_to_calib_format as movieToCalFmt


def test_create_cal_1():
    x_size = 110
    y_size = 100
    gain = 2.0
    offset = 100.0
    read_noise = 1.0
    [cam_offset, cam_variance, cam_gain, cam_rqe] = camera.createSCMOSCalibration(x_size,
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
    assert numpy.allclose(cam_rqe, numpy.ones((x_size, y_size)))

    
def test_create_cal_2():
    x_size = 110
    y_size = 100
    gain = 2.0
    offset = 100.0
    read_noise = 0.0
    [cam_offset, cam_variance, cam_gain, cam_rqe] = camera.createSCMOSCalibration(x_size,
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
    assert numpy.allclose(cam_rqe, numpy.ones((x_size, y_size)))
    
    
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


def test_cam_cal_1():
    """
    Calibration file format 0.
    """
    size = (12,10)
    cam_gain = 1.5 * numpy.ones(size)
    cam_offset = 1000.0 * numpy.ones(size)
    cam_var = 2.0 * numpy.ones(size)
    n_frames = 20000
    
    # Create calibration files.
    scmos_files = []
    for i, name in enumerate(["dark.npy", "light1.npy", "light2.npy", "light3.npy", "light4.npy"]):
        f_name = storm_analysis.getPathOutputTest(name)
        scmos_files.append(f_name)
        
        mean = i * 500 * cam_gain
        var = mean * cam_gain + cam_var
        mean += cam_offset

        N = mean * n_frames
        NN = (var + mean*mean) * n_frames
        
        numpy.save(f_name, [numpy.array([n_frames]), N, NN])

    # Check.
    [cal_offset, cal_var, cal_gain, cal_rqe] = camCal.cameraCalibration(scmos_files,
                                                                        show_fit_plots = False,
                                                                        show_mean_plots = False)

    assert(numpy.allclose(cal_offset, cam_offset))
    assert(numpy.allclose(cal_var, cam_var))
    assert(numpy.allclose(cal_gain, cam_gain))
    assert(numpy.allclose(cal_rqe, numpy.ones(size)))


def test_cam_cal_2():
    """
    Calibration file format 1.
    """
    size = (12,10)
    cam_gain = 1.5 * numpy.ones(size)
    cam_offset = 1000.0 * numpy.ones(size)
    cam_var = 2.0 * numpy.ones(size)
    n_frames = 20000
    
    # Create calibration files.
    scmos_files = []
    for i, name in enumerate(["dark.npy", "light1.npy", "light2.npy", "light3.npy", "light4.npy"]):
        f_name = storm_analysis.getPathOutputTest(name)
        scmos_files.append(f_name)
        
        mean = i * 500 * cam_gain
        var = mean * cam_gain + cam_var
        mean += cam_offset

        N = mean * n_frames
        NN = (var + mean*mean) * n_frames

        mean_mean = numpy.zeros(n_frames) + numpy.mean(mean)
        numpy.save(f_name, [mean_mean, N, NN])

    # Check.
    [cal_offset, cal_var, cal_gain, cal_rqe] = camCal.cameraCalibration(scmos_files,
                                                                        show_fit_plots = False,
                                                                        show_mean_plots = False)

    assert(numpy.allclose(cal_offset, cam_offset))
    assert(numpy.allclose(cal_var, cam_var))
    assert(numpy.allclose(cal_gain, cam_gain))
    assert(numpy.allclose(cal_rqe, numpy.ones(size)))


def test_cam_cal_3():
    """
    Calibration file format 2.
    """
    size = (12,10)
    cam_gain = 1.5 * numpy.ones(size)
    cam_offset = 1000.0 * numpy.ones(size)
    cam_var = 2.0 * numpy.ones(size)
    n_frames = 20000
    
    # Create calibration files.
    scmos_files = []
    for i, name in enumerate(["dark.npy", "light1.npy", "light2.npy", "light3.npy", "light4.npy"]):
        f_name = storm_analysis.getPathOutputTest(name)
        scmos_files.append(f_name)
        
        mean = i * 500 * cam_gain
        var = mean * cam_gain + cam_var
        mean += cam_offset

        N = mean * n_frames
        NN = (var + mean*mean) * n_frames

        mean_mean = numpy.zeros(n_frames) + numpy.mean(mean)
        roi_dict = {"x_start" : 1, "y_start" : 2}
        numpy.save(f_name, [mean_mean, N, NN])

    # Check.
    [cal_offset, cal_var, cal_gain, cal_rqe] = camCal.cameraCalibration(scmos_files,
                                                                        show_fit_plots = False,
                                                                        show_mean_plots = False)

    assert(numpy.allclose(cal_offset, cam_offset))
    assert(numpy.allclose(cal_var, cam_var))
    assert(numpy.allclose(cal_gain, cam_gain))
    assert(numpy.allclose(cal_rqe, numpy.ones(size)))


def test_bad_pixel():
    """
    Test bad pixel handling.
    """
    size = (12,10)
    cam_gain = 1.5 * numpy.ones(size)
    cam_offset = 1000.0 * numpy.ones(size)
    cam_var = 2.0 * numpy.ones(size)
    n_frames = 20000
    
    # Create calibration files.
    scmos_files = []
    for i, name in enumerate(["dark.npy", "light1.npy", "light2.npy", "light3.npy", "light4.npy"]):
        f_name = storm_analysis.getPathOutputTest(name)
        scmos_files.append(f_name)
        
        mean = i * 500 * cam_gain
        mean[4,5] = 0.0
        
        var = mean * cam_gain + cam_var
        mean += cam_offset

        N = mean * n_frames
        NN = (var + mean*mean) * n_frames

        mean_mean = numpy.zeros(n_frames) + numpy.mean(mean)
        roi_dict = {"x_start" : 1, "y_start" : 2}
        numpy.save(f_name, [mean_mean, N, NN])

    # Check.
    [cal_offset, cal_var, cal_gain, cal_rqe] = camCal.cameraCalibration(scmos_files,
                                                                        show_fit_plots = False,
                                                                        show_mean_plots = False)

    cam_gain[4,5] = 1.0
    
    assert(numpy.allclose(cal_offset, cam_offset))
    assert(numpy.allclose(cal_var, cam_var))
    assert(numpy.allclose(cal_gain, cam_gain))
    
    
if (__name__ == "__main__"):
    test_create_cal_1()
    test_create_cal_2()
    test_mtcf_1()
    test_cam_cal_1()
    test_cam_cal_2()
    test_cam_cal_3()
    test_bad_pixel()
    
    
