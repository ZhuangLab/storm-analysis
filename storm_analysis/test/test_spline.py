#!/usr/bin/env python
"""
Tests of splines.
"""
import numpy
import pickle
import random
import sys

import storm_analysis

import storm_analysis.spliner.cubic_spline_c as cubicSplineC
import storm_analysis.spliner.spline2D as spline2D
import storm_analysis.spliner.spline3D as spline3D


reps = 1000

def test_psf_2D_f():

    # Only test for Python3 due to pickle incompatibility issues which I am tired
    # of trying to deal with.
    #
    if (sys.version_info < (3, 0)):
        return
    
    spline_filename = storm_analysis.getData("test/data/test_spliner_psf_2d.spline")
    with open(spline_filename, "rb") as fp:
        spline_data = pickle.load(fp)

    py_spline = spline2D.Spline2D(spline_data["spline"], spline_data["coeff"])
    c_spline = cubicSplineC.CSpline2D(py_spline)

    size = py_spline.getSize() - 1.0e-6

    for i in range(reps):
        x = random.uniform(1.0e-6, size)
        y = random.uniform(1.0e-6, size)
        #print("{0:.3f} {1:.3f}".format(py_spline.f(x, y), c_spline.f(x, y)))
        assert (abs(py_spline.f(x, y) - c_spline.f(x, y)) < 1.0e-6)

def test_psf_2D_dx():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_filename = storm_analysis.getData("test/data/test_spliner_psf_2d.spline")
    with open(spline_filename, "rb") as fp:
        spline_data = pickle.load(fp)

    py_spline = spline2D.Spline2D(spline_data["spline"], spline_data["coeff"])
    c_spline = cubicSplineC.CSpline2D(py_spline)

    size = py_spline.getSize() - 1.0e-6

    for i in range(reps):
        x = random.uniform(1.0e-6, size)
        y = random.uniform(1.0e-6, size)
        #print("{0:.3f} {1:.3f}".format(py_spline.dxf(x, y), c_spline.dxf(x, y)))
        assert (abs(py_spline.dxf(x, y) - c_spline.dxf(x, y)) < 1.0e-6)        

def test_psf_2D_dy():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_filename = storm_analysis.getData("test/data/test_spliner_psf_2d.spline")
    with open(spline_filename, "rb") as fp:
        spline_data = pickle.load(fp)

    py_spline = spline2D.Spline2D(spline_data["spline"], spline_data["coeff"])
    c_spline = cubicSplineC.CSpline2D(py_spline)

    size = py_spline.getSize() - 1.0e-6

    for i in range(reps):
        x = random.uniform(1.0e-6, size)
        y = random.uniform(1.0e-6, size)
        #print("{0:.3f} {1:.3f}".format(py_spline.dyf(x, y), c_spline.dyf(x, y)))
        assert (abs(py_spline.dyf(x, y) - c_spline.dyf(x, y)) < 1.0e-6)


def test_psf_3D_f():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_filename = storm_analysis.getData("test/data/test_spliner_psf.spline")
    with open(spline_filename, "rb") as fp:
        spline_data = pickle.load(fp)

    py_spline = spline3D.Spline3D(spline_data["spline"], spline_data["coeff"])
    c_spline = cubicSplineC.CSpline3D(py_spline)

    size = py_spline.getSize() - 1.0e-6

    for i in range(reps):
        x = random.uniform(1.0e-6, size)
        y = random.uniform(1.0e-6, size)
        z = random.uniform(1.0e-6, size)
        #print("{0:.3f} {1:.3f}".format(py_spline.f(x, y), c_spline.f(x, y)))
        assert (abs(py_spline.f(x, y, z) - c_spline.f(x, y, z)) < 1.0e-6)

def test_psf_3D_dx():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_filename = storm_analysis.getData("test/data/test_spliner_psf.spline")
    with open(spline_filename, "rb") as fp:
        spline_data = pickle.load(fp)

    py_spline = spline3D.Spline3D(spline_data["spline"], spline_data["coeff"])
    c_spline = cubicSplineC.CSpline3D(py_spline)

    size = py_spline.getSize() - 1.0e-6

    for i in range(reps):
        x = random.uniform(1.0e-6, size)
        y = random.uniform(1.0e-6, size)
        z = random.uniform(1.0e-6, size)
        #print("{0:.3f} {1:.3f}".format(py_spline.dxf(x, y), c_spline.dxf(x, y)))
        assert (abs(py_spline.dxf(x, y, z) - c_spline.dxf(x, y, z)) < 1.0e-6)

def test_psf_3D_dy():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_filename = storm_analysis.getData("test/data/test_spliner_psf.spline")
    with open(spline_filename, "rb") as fp:
        spline_data = pickle.load(fp)

    py_spline = spline3D.Spline3D(spline_data["spline"], spline_data["coeff"])
    c_spline = cubicSplineC.CSpline3D(py_spline)

    size = py_spline.getSize() - 1.0e-6

    for i in range(reps):
        x = random.uniform(1.0e-6, size)
        y = random.uniform(1.0e-6, size)
        z = random.uniform(1.0e-6, size)
        #print("{0:.3f} {1:.3f}".format(py_spline.dyf(x, y), c_spline.dyf(x, y)))
        assert (abs(py_spline.dyf(x, y, z) - c_spline.dyf(x, y, z)) < 1.0e-6)

def test_psf_3D_dz():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_filename = storm_analysis.getData("test/data/test_spliner_psf.spline")
    with open(spline_filename, "rb") as fp:
        spline_data = pickle.load(fp)

    py_spline = spline3D.Spline3D(spline_data["spline"], spline_data["coeff"])
    c_spline = cubicSplineC.CSpline3D(py_spline)

    size = py_spline.getSize() - 1.0e-6

    for i in range(reps):
        x = random.uniform(1.0e-6, size)
        y = random.uniform(1.0e-6, size)
        z = random.uniform(1.0e-6, size)
        #print("{0:.3f} {1:.3f}".format(py_spline.dzf(x, y), c_spline.dzf(x, y)))
        assert (abs(py_spline.dzf(x, y, z) - c_spline.dzf(x, y, z)) < 1.0e-6)


if (__name__ == "__main__"):
    test_psf_2D_f()
    test_psf_2D_dx()
    test_psf_2D_dy()
    test_psf_3D_f()
    test_psf_3D_dx()
    test_psf_3D_dy()
    test_psf_3D_dz()
