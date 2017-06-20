#!/usr/bin/env python
"""
Tests for sa_library.matched_filter_c.py

Hazen 06/17
"""
import numpy

import storm_analysis.sa_library.matched_filter_c as matchedFilterC
import storm_analysis.simulator.draw_gaussians_c as dg


def test_matched_filter1():
    x_size = 80
    y_size = 90

    objects = numpy.zeros((1, 5))

    # Make filter with unit sum.    
    objects[0,:] = [x_size/2, y_size/2, 1.0, 1.0, 1.0]
    psf = dg.drawGaussians((x_size, y_size), objects)
    psf = psf/numpy.sum(psf)
    flt = matchedFilterC.MatchedFilter(psf)

    image = numpy.zeros((x_size, y_size))
    image[int(x_size/2), int(y_size/2)] = 2.0

    conv = flt.convolve(image)
    assert(abs(numpy.sum(image) - numpy.sum(conv)) < 1.0e-6)


if (__name__ == "__main__"):
    test_matched_filter1()

