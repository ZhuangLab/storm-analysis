#!/usr/bin/env python
"""
Tests for sa_library.imagecorrelation
"""
import numpy
import matplotlib
import matplotlib.pyplot as pyplot
import tifffile

import storm_analysis.sa_library.imagecorrelation as imgCorr
import storm_analysis.simulator.draw_gaussians_c as dg


def test_align3d_1():
    """
    Test on uniform images.
    """
    image = numpy.ones((10,12,14))
    a3d = imgCorr.Align3DProduct(image, xy_margin = 1, z_margin = 1)
    a3d.setOtherImage(image)

    # Check for the same value regardless of displacement.
    val = a3d.func(numpy.zeros(3))
    check = numpy.array([[0.0, 0.0, 0.1],
                         [0.0, 0.1, 0.0],
                         [0.0, 0.1, 0.1],
                         [0.1, 0.0, 0.1],
                         [0.1, -0.1, 0.0],
                         [0.1, -0.1, 0.1]])
    for i in range(check.shape[0]):
        assert(abs(val - a3d.func(check[i,:])) < 1.0e-6)

    # Check that derivatives are zero.
    for elt in [-0.1, 0.0, 0.1]:
        assert(abs(a3d.jacobian(numpy.array([elt, 0.0, 0.0]))[0]) < 1.0e-6)
        assert(abs(a3d.jacobian(numpy.array([0.0, elt, 0.0]))[1]) < 1.0e-6)
        assert(abs(a3d.jacobian(numpy.array([0.0, 0.0, elt]))[2]) < 1.0e-6)


def test_align3d_2():
    """
    Test displacement
    """
    image1 = dg.drawGaussiansXYZ((12,13,14),
                                 numpy.array([6.0]),
                                 numpy.array([6.5]),
                                 numpy.array([7.0]))
    a3d = imgCorr.Align3DProduct(image1, xy_margin = 1, z_margin = 1)
    a3d.setOtherImage(image1)
    val = a3d.func(numpy.zeros(3))
    
    image2 = dg.drawGaussiansXYZ((12,13,14),
                                 numpy.array([6.1]),
                                 numpy.array([6.6]),
                                 numpy.array([6.9]))    
    a3d.setOtherImage(image2)

    
    assert(a3d.func(numpy.array([0.0, 0.0, 0.0])) < (val - 0.01))
    assert(a3d.func(numpy.array([0.1, 0.1, -0.1])) < (val - 0.01))
    assert(abs(a3d.func(numpy.array([-0.1, -0.1, 0.1])) - val) < 1.0e-4)


def test_align3d_3():
    """
    Test on Gaussians.
    """
    image = dg.drawGaussiansXYZ((12,13,14),
                                numpy.array([6.0]),
                                numpy.array([6.5]),
                                numpy.array([7.0]))
    
    a3d = imgCorr.Align3DProduct(image, xy_margin = 1, z_margin = 1)
    a3d.setOtherImage(image)

    # Check X derivatives.
    for elt in [-0.1, 0.0, 0.1]:
        dx = 1.0e-6
        f1 = a3d.func(numpy.array([elt + dx, 0.0, 0.0]))
        f2 = a3d.func(numpy.array([elt, 0.0, 0.0]))
        est = (f1-f2)/dx
        exact = a3d.jacobian(numpy.array([elt, 0.0, 0.0]))[0]
        assert(abs(est - exact) < 1.0e-3)
    
    # Check Y derivatives.
    for elt in [-0.1, 0.0, 0.1]:
        dx = 1.0e-6
        f1 = a3d.func(numpy.array([0.0, elt + dx, 0.0]))
        f2 = a3d.func(numpy.array([0.0, elt, 0.0]))
        est = (f1-f2)/dx
        exact = a3d.jacobian(numpy.array([0.0, elt, 0.0]))[1]
        assert(abs(est - exact) < 1.0e-3)

    # Check Z derivatives.
    for elt in [-0.1, 0.0, 0.1]:
        dx = 1.0e-6
        f1 = a3d.func(numpy.array([0.0, 0.0, elt + dx]))
        f2 = a3d.func(numpy.array([0.0, 0.0, elt]))
        est = (f1-f2)/dx
        exact = a3d.jacobian(numpy.array([0.0, 0.0, elt]))[2]
        assert(abs(est - exact) < 1.0e-3)


def test_align3d_4():
    """
    Test alignment.
    """
    dx = numpy.array([0.1, 0.15, -0.08])
    #dx = numpy.array([0.1, 0.01, 0.001])
    image1 = dg.drawGaussiansXYZ((12,13,14),
                                 numpy.array([6.0]),
                                 numpy.array([6.5]),
                                 numpy.array([7.0]))
#    a3d = imgCorr.Align3DProductLM(image1, xy_margin = 1, z_margin = 1)
    a3d = imgCorr.Align3DProductNewtonCG(image1, xy_margin = 1, z_margin = 1)
    
    image2 = dg.drawGaussiansXYZ((12,13,14),
                                 numpy.array([6.0 - dx[0]]),
                                 numpy.array([6.5 - dx[1]]),
                                 numpy.array([7.0 - dx[2]]))    
    a3d.setOtherImage(image2)

    # Check that we find the proper alignment.
    [disp, success, fun, status] = a3d.maximize()
    assert(success)
    assert(numpy.allclose(disp, dx, rtol = 0.01, atol = 0.01))

    # Check returning an aligned other.
    [aligned, q_score, disp] = a3d.align()
    assert(q_score > 10.0)
    assert(numpy.allclose(image1, aligned, rtol = 1.0e-2, atol = 1.0e-2))
    
    
if (__name__ == "__main__"):
    test_align3d_1()
    test_align3d_2()
    test_align3d_3()
    test_align3d_4()
