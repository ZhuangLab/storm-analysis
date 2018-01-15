#!/usr/bin/env python
"""
Test HDF5 to image conversion.
"""
import h5py
import numpy
import os

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


if (__name__ == "__main__"):
    test_hdf5_to_image_1()
