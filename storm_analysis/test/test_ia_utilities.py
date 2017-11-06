#!/usr/bin/env python
"""
Tests some aspect of ia_utilities.
"""
import numpy
import tifffile

import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.simulator.draw_gaussians_c as dg


def test_ia_util_1():
    """
    Test finding peaks in an empty image.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(images)
    assert(x.size == 0)

def test_ia_util_2():
    """
    Test finding peaks in an image.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    # Above threshold, unequal height.
    images[0][10,11] = 1.1
    images[0][10,12] = 1.5

    # Abover threshold, unequal height.
    images[0][15,8] = 1.5
    images[0][15,9] = 1.5

    # Below threshold.
    images[0][20,11] = 0.5
    
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(images)
    print(x, y, z)
    
if (__name__ == "__main__"):
    test_ia_util_1()
    test_ia_util_2()


