#!/usr/bin/env python
"""
Tests for spliner.measure_psf_utils
"""
import numpy
import random
import scipy
import tifffile

import storm_analysis

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.parameters as params
import storm_analysis.simulator.draw_gaussians_c as dg
import storm_analysis.spliner.measure_psf_utils as mPSFUtils


def test_mzia_1():
    """
    Test z_index array creation.
    """
    z_offsets = numpy.array([[1, -0.2501],
                             [1, -0.21],
                             [1, -0.151],
                             [1, -0.149],
                             [1, 0.0],
                             [1, 2.49],
                             [1, 2.5],
                             [0, 0.0]])
    z_index = mPSFUtils.makeZIndexArray(z_offsets, 0.2, 0.1)
    print(z_index)
    assert(numpy.allclose(z_index, numpy.array([-1, 0, 0, 1, 2, -1, -1, -1])))


def test_mspb_1():
    """
    Test (single) PSF measurement, no drift.
    """
    
    # Make test movie.
    x = 7.2
    y = 9.8
    psf_movie = storm_analysis.getPathOutputTest("psf_movie.tif")
    image = 1000.0*dg.drawGaussiansXY((20, 20), numpy.array([x]), numpy.array([y]))
    with tifffile.TiffWriter(psf_movie) as tf:
        for i in range(6):
            tf.save(image.astype(numpy.float32))

    # Parameters.
    p = params.ParametersDAO()
    p.changeAttr("camera_gain", 1.0)
    p.changeAttr("camera_offset", 0.0)

    # Frame reader.
    frdr = analysisIO.FrameReaderStd(movie_file = psf_movie,
                                     parameters = p)

    z_index = numpy.array([0, 1, 2, 2, -1, -1])
    [psf, samples] = mPSFUtils.measureSinglePSFBeads(frdr, z_index, 6, x, y, zoom = 2)
    
    assert(numpy.allclose(samples, numpy.array([1,1,2])))
    for i in range(1,psf.shape[0]):
        assert(numpy.allclose(psf[0,:,:],psf[i,:,:]/samples[i]))

    if False:
        with tifffile.TiffWriter("psf.tif") as tf:
            for i in range(psf.shape[0]):
                tf.save(psf[i,:,:].astype(numpy.float32))


def test_mspb_2():
    """
    Test (single) PSF measurement, no drift, recentering.

    The maximum relative difference is typically on the order of 2%.
    """
    
    # Make test movie.
    im_max = 1000.0
    x = 7.0 + numpy.random.uniform(size = 10)
    y = 11.0 + numpy.random.uniform(size = 10)
    
    psf_movie = storm_analysis.getPathOutputTest("psf_movie.tif")
    with tifffile.TiffWriter(psf_movie) as tf:
        for i in range(x.size):
            image = dg.drawGaussiansXY((20, 20),
                                       numpy.array([x[i]]),
                                       numpy.array([y[i]]))
            image = image * im_max
            tf.save(image.astype(numpy.float32))

    # Parameters.
    p = params.ParametersDAO()
    p.changeAttr("camera_gain", 1.0)
    p.changeAttr("camera_offset", 0.0)

    # Frame reader.
    frdr = analysisIO.FrameReaderStd(movie_file = psf_movie,
                                     parameters = p)
    
    z_index = numpy.zeros(x.size).astype(numpy.int) - 1
    z_index[0] = 0
    [psf0, samples] = mPSFUtils.measureSinglePSFBeads(frdr, z_index, 6, x[0], y[0], zoom = 2)

    for i in range(1,x.size):
        z_index = numpy.zeros(x.size).astype(numpy.int) - 1
        z_index[i] = 0
        [psf, samples] = mPSFUtils.measureSinglePSFBeads(frdr, z_index, 6, x[i], y[i], zoom = 2)
        assert(numpy.max(numpy.abs(psf0 - psf)/numpy.max(psf)) < 0.05)


def test_mspb_3():
    """
    Test (single) PSF measurement with drift.
    """
    
    # Make test movie.
    im_max = 1000.0
    n_pts = 10
    x = 7.0
    y = 11.0
    drift_xy = numpy.random.uniform(size = (n_pts, 2))
    
    psf_movie = storm_analysis.getPathOutputTest("psf_movie.tif")
    with tifffile.TiffWriter(psf_movie) as tf:
        for i in range(n_pts):
            image = dg.drawGaussiansXY((20, 20),
                                       numpy.array([x + drift_xy[i][0]]),
                                       numpy.array([y + drift_xy[i][1]]))
            image = image * im_max
            tf.save(image.astype(numpy.float32))

    # Parameters.
    p = params.ParametersDAO()
    p.changeAttr("camera_gain", 1.0)
    p.changeAttr("camera_offset", 0.0)

    # Frame reader.
    frdr = analysisIO.FrameReaderStd(movie_file = psf_movie,
                                     parameters = p)
    
    z_index = numpy.zeros(n_pts).astype(numpy.int) - 1
    z_index[0] = 0
    [psf0, samples] = mPSFUtils.measureSinglePSFBeads(frdr, z_index, 6, x + drift_xy[0][0], y + drift_xy[0][1], zoom = 2)

    for i in range(1, n_pts):
        z_index = numpy.zeros(n_pts).astype(numpy.int) - 1
        z_index[i] = 0
        [psf, samples] = mPSFUtils.measureSinglePSFBeads(frdr, z_index, 6, x, y, drift_xy = drift_xy, zoom = 2)
        assert(numpy.max(numpy.abs(psf0 - psf)/numpy.max(psf)) < 0.05)


def test_align_psfs_1():
    """
    Test alignment of multiple PSFs to each other.
    """
    pos = [5.0, 6.0, 7.0]
    maxd = 3.0

    psfs = []
    for i in range(7):
        psfs.append(dg.drawGaussiansXYZ((12,13,14),
                                        numpy.array([pos[0] + int(i/maxd)]),
                                        numpy.array([pos[1] + int(i/maxd)]),
                                        numpy.array([pos[2] + int(i/maxd)])))
    
    [average_psf, i_score] = mPSFUtils.alignPSFs(psfs)
    assert(i_score > 2.1)

    if False:
        with tifffile.TiffWriter("psf.tif") as tf:
            tf.save(mPSFUtils.averagePSF(psfs)[6,:,:].astype(numpy.float32))
            tf.save(average_psf[6,:,:].astype(numpy.float32))


def test_zscaler_1():
    """
    Test ZScaler for floating point round off issues.
    """
    zs = mPSFUtils.ZScaler(0.6, 0.2)
    assert(zs.getMaxZ() == 7)

    diff = 1.0e-24
    zs = mPSFUtils.ZScaler(0.6 + diff, 0.2 - diff)
    assert(zs.getMaxZ() == 7)

    zs = mPSFUtils.ZScaler(0.6 - diff, 0.2 + diff)
    assert(zs.getMaxZ() == 7)


def test_extract_roi_1():
    """
    Test range checking in extractROI.
    """
    frame = numpy.zeros((100,100))

    aoi_size = 10

    # X test 1.
    im = mPSFUtils.extractAOI(frame, aoi_size, 10, 20)
    assert(im.shape[0] == 2*aoi_size)
    assert(im.shape[1] == 2*aoi_size)

    okay = False
    try:
        mPSFUtils.extractAOI(frame, aoi_size, 9, 20)
    except AssertionError:
        okay = True
    assert okay

    # Y test 1.
    im = mPSFUtils.extractAOI(frame, aoi_size, 20, 10)
    assert(im.shape[0] == 2*aoi_size)
    assert(im.shape[1] == 2*aoi_size)

    okay = False
    try:
        mPSFUtils.extractAOI(frame, aoi_size, 20, 9)
    except AssertionError:
        okay = True
    assert okay

    # X test 2.
    im = mPSFUtils.extractAOI(frame, aoi_size, 90, 20)
    assert(im.shape[0] == 2*aoi_size)
    assert(im.shape[1] == 2*aoi_size)

    okay = False
    try:
        mPSFUtils.extractAOI(frame, aoi_size, 91, 20)
    except AssertionError:
        okay = True
    assert okay

    # Y test 2.
    im = mPSFUtils.extractAOI(frame, aoi_size, 20, 90)
    assert(im.shape[0] == 2*aoi_size)
    assert(im.shape[1] == 2*aoi_size)

    okay = False
    try:
        mPSFUtils.extractAOI(frame, aoi_size, 20, 91)
    except AssertionError:
        okay = True
    assert okay


if (__name__ == "__main__"):
    test_mzia_1()
    test_mspb_1()
    test_mspb_2()
    test_mspb_3()
    test_align_psfs_1()
    test_zscaler_1()
    
