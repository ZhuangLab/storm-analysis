#!/usr/bin/env python
"""
Tests for sa_library.matched_filter_c.py

Hazen 06/17
"""
import numpy

import storm_analysis.sa_library.matched_filter_c as matchedFilterC
import storm_analysis.simulator.draw_gaussians_c as dg


def test_matched_filter1():
    """
    Verify that the filter results are normalized correctly.
    """
    x_size = 80
    y_size = 90

    objects = numpy.zeros((1, 5))

    # Make filter with unit sum.
    objects[0,:] = [x_size/2, y_size/2, 1.0, 1.0, 1.0]
    psf = dg.drawGaussians((x_size, y_size), objects)
    psf = psf/numpy.sum(psf)
    flt = matchedFilterC.MatchedFilter(psf)

    for i in range(1,5):
        image = numpy.zeros((x_size, y_size))
        image[int(x_size/2), int(y_size/2)] = float(i)
        conv = flt.convolve(image)
        assert(abs(numpy.sum(image) - numpy.sum(conv)) < 1.0e-6)

def test_matched_filter2():
    """
    Test recovering the original height from the convolved image.
    """
    x_size = 80
    y_size = 90

    objects = numpy.zeros((1, 5))

    # Make filter with unit sum.    
    objects[0,:] = [x_size/2, y_size/2, 1.0, 2.0, 2.0]
    psf = dg.drawGaussians((x_size, y_size), objects)
    psf_norm = psf/numpy.sum(psf)
    flt = matchedFilterC.MatchedFilter(psf_norm)

    rescale = 1.0/numpy.sum(psf * psf_norm)

    # Create object with height 10 and the same shape as the filter.
    height = 10.0
    objects[0,:] = [x_size/2, y_size/2, height, 2.0, 2.0]
    image = dg.drawGaussians((x_size, y_size), objects)

    # Apply filter.
    conv = flt.convolve(image)

    # Verify that final height is 'close enough'.
    assert (abs(numpy.amax(conv) * rescale - height)/height < 1.0e-2)
    

if (__name__ == "__main__"):
    test_matched_filter1()
    test_matched_filter2()

