#!/usr/bin/env python
"""
Test HDF5 to image conversion.
"""
import h5py
import numpy
import os
import tifffile


import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_library.writeinsight3 as i3w
import storm_analysis.sa_utilities.hdf5_to_image as hdf5ToImage


def test_hdf5_to_image_1():
    """
    Test simple 2D image (using localizations).
    """
    peaks = {"x" : numpy.array([10.0]),
             "y" : numpy.array([20.0])}

    h5_name = storm_analysis.getPathOutputTest("test_hdf5_to_image.hdf5")
    storm_analysis.removeFile(h5_name)

    # Write data (as peaks).
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(40,30,1,"")
        h5.addLocalizations(peaks, 0)

    # Render image.
    image = hdf5ToImage.render2DImage(h5_name, scale = 1)

    assert(image.shape[0] == 30)
    assert(image.shape[1] == 40)
    assert(image[10,20] == 0)
    assert(image[20,10] == 1)


def test_hdf5_to_image_2():
    """
    Test simple 2D image (using tracks).
    """
    tracks = {"x" : numpy.array([10.0]),
              "y" : numpy.array([20.0])}

    h5_name = storm_analysis.getPathOutputTest("test_hdf5_to_image.hdf5")
    storm_analysis.removeFile(h5_name)

    # Write data (as peaks).
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(40,30,1,"")
        h5.addTracks(tracks)

    # Render image.
    image = hdf5ToImage.render2DImage(h5_name, scale = 1)

    assert(image.shape[0] == 30)
    assert(image.shape[1] == 40)
    assert(image[10,20] == 0)
    assert(image[20,10] == 1)


def test_hdf5_to_image_3():
    """
    Test simple 2D image versus gaussian image.
    """
    tracks = {"x" : numpy.array([10]),
              "y" : numpy.array([20])}

    h5_name = storm_analysis.getPathOutputTest("test_hdf5_to_image.hdf5")
    storm_analysis.removeFile(h5_name)

    # Write data (as peaks).
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(40,30,1,"")
        h5.addTracks(tracks)

    # Render image.
    im_grid = hdf5ToImage.render2DImage(h5_name, scale = 1)
    im_gauss = hdf5ToImage.render2DImage(h5_name, scale = 1, sigma = 0.25)

    assert(numpy.allclose(im_grid, im_gauss, atol = 1.0e-3, rtol = 1.0e-3))
    

def test_hdf5_to_image_4():
    """
    Test category (using tracks).
    """
    tracks = {"category" : [1, 2],
              "x" : numpy.array([10.0, 20.0]),
              "y" : numpy.array([20.0, 10.0])}

    h5_name = storm_analysis.getPathOutputTest("test_hdf5_to_image.hdf5")
    storm_analysis.removeFile(h5_name)

    # Write data (as peaks).
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(40,30,1,"")
        h5.addTracks(tracks)

    # Render image.
    image = hdf5ToImage.render2DImage(h5_name, category = 1, scale = 1)

    assert(image.shape[0] == 30)
    assert(image.shape[1] == 40)
    assert(image[10,20] == 0)
    assert(image[20,10] == 1)


def test_hdf5_to_image_5():
    """
    Test 3D rendering.
    """
    tracks = {"x" : numpy.array([10.0, 10.0, 10.0]),
              "y" : numpy.array([20.0, 20.0, 20.0]),
              "z" : numpy.array([-0.2, 0.0, 0.2])}

    h5_name = storm_analysis.getPathOutputTest("test_hdf5_to_image.hdf5")
    storm_analysis.removeFile(h5_name)

    # Write data (as peaks).
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(40,30,1,"")
        h5.addTracks(tracks)

    # Render image.
    images = hdf5ToImage.render3DImage(h5_name, [-0.3, -0.1, 0.1, 0.3], scale = 1)

    assert(images[0].shape[0] == 30)
    assert(images[0].shape[1] == 40)
    assert(numpy.allclose(images[0], images[1]))
    assert(numpy.allclose(images[0], images[2]))

    
if (__name__ == "__main__"):
    test_hdf5_to_image_1()
    test_hdf5_to_image_2()
    test_hdf5_to_image_3()
    test_hdf5_to_image_4()
    test_hdf5_to_image_5()
