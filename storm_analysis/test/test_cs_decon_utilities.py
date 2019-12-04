#!/usr/bin/env python
"""
Tests some aspect of cs_decon_utilities.
"""
import numpy

import storm_analysis.sa_library.cs_decon_utilities_c as csDeconUtilsC


def test_cs_utils_1():
    
    # Create test image with some non-zero data.
    test_image = numpy.zeros((10,10,3))

    test_image[5,0,0] = 1.0
    
    test_image[1,1,0] = 1.0
    
    test_image[5,5,0] = 1.0
    test_image[5,5,1] = 2.0

    test_image[7,6,1] = 1.0
    test_image[7,7,1] = 2.0
    test_image[7,8,1] = 1.0

    [labels, counts] = csDeconUtilsC.label(test_image, 0.1, 0)
    assert(counts == 3)

    peaks = csDeconUtilsC.getPeaks(test_image, 0.1, 0)
    assert(numpy.allclose(peaks[:,1], numpy.array([1, 5, 7])))
    assert(numpy.allclose(peaks[:,2], numpy.array([1, 5, 7])))


if (__name__ == "__main__"):
    test_cs_utils_1()
    
