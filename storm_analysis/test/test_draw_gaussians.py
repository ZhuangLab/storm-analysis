#!/usr/bin/env python
"""
Tests of simulator.draw_gaussians_c
"""
import numpy
import tifffile

import storm_analysis.simulator.draw_gaussians_c as dg


def test_dgxy_1():
    image = dg.drawGaussiansXY((20,20),
                               numpy.array([10.0]),
                               numpy.array([10.0]))

    dx = numpy.arange(-10.0,10,1.0)
    g_slice = numpy.exp(-dx*dx/2.0)

    assert(numpy.allclose(g_slice, image[10,:], atol = 1.0e-4))
    assert(numpy.allclose(g_slice, image[:,10], atol = 1.0e-4))


def test_dgxy_2():
    image = dg.drawGaussiansXYZ((20,20,20),
                                numpy.array([10.0]),
                                numpy.array([10.0]),
                                numpy.array([10.0]))

    dx = numpy.arange(-10.0,10,1.0)
    g_slice = numpy.exp(-dx*dx/2.0)

    assert(numpy.allclose(g_slice, image[10,10,:], atol = 1.0e-4))
    assert(numpy.allclose(g_slice, image[10,:,10], atol = 1.0e-4))
    assert(numpy.allclose(g_slice, image[:,10,10], atol = 1.0e-4))
    
    
if (__name__ == "__main__"):
    test_dgxy_1()
    test_dgxy_2()

    
