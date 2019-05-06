#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.simulator.pupil_math as pupilMath


def test_pupil_math_1():
    """
    Test GeometryC, intensity, no scaling.
    """
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    geo_c = pupilMath.GeometryC(20, 0.1, 0.6, 1.5, 1.4)

    pf = geo.createFromZernike(1.0,  [[1.3, -1, 3], [1.3, -2, 2]])    
    z_vals = numpy.linspace(-1.0,1.0,10)

    psf_py = geo.pfToPSF(pf, z_vals)
    psf_c = geo_c.pfToPSF(pf, z_vals)
    
    assert numpy.allclose(psf_c, psf_py)

def test_pupil_math_2():
    """
    Test GeometryC, complex values, no scaling.
    """
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    geo_c = pupilMath.GeometryC(20, 0.1, 0.6, 1.5, 1.4)

    pf = geo.createFromZernike(1.0,  [[1.3, -1, 3], [1.3, -2, 2]])    
    z_vals = numpy.linspace(-1.0,1.0,10)

    psf_py = geo.pfToPSF(pf, z_vals, want_intensity = False)
    psf_c = geo_c.pfToPSF(pf, z_vals, want_intensity = False)
    
    assert numpy.allclose(psf_c, psf_py)

def test_pupil_math_3():
    """
    Test GeometryC, intensity, scaling.
    """
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    geo_c = pupilMath.GeometryC(20, 0.1, 0.6, 1.5, 1.4)

    pf = geo.createFromZernike(1.0,  [[1.3, -1, 3], [1.3, -2, 2]])    
    z_vals = numpy.linspace(-1.0,1.0,10)

    gsf = geo.gaussianScalingFactor(1.8)
    psf_py = geo.pfToPSF(pf, z_vals, scaling_factor = gsf)
    psf_c = geo_c.pfToPSF(pf, z_vals, scaling_factor = gsf)
    
    assert numpy.allclose(psf_c, psf_py)

def test_pupil_math_4():
    """
    Test GeometryCVectorial, intensity, no scaling.
    """
    geo = pupilMath.GeometryVectorial(20, 0.1, 0.6, 1.5, 1.4)
    geo_c = pupilMath.GeometryCVectorial(20, 0.1, 0.6, 1.5, 1.4)

    pf = geo.createFromZernike(1.0,  [[1.3, -1, 3], [1.3, -2, 2]])    
    z_vals = numpy.linspace(-1.0,1.0,10)

    psf_py = geo.pfToPSF(pf, z_vals)
    psf_c = geo_c.pfToPSF(pf, z_vals)
    
    assert numpy.allclose(psf_c, psf_py)

def test_pupil_math_5():
    """
    Test GeometryCVectorial, intensity, scaling.
    """
    geo = pupilMath.GeometryVectorial(20, 0.1, 0.6, 1.5, 1.4)
    geo_c = pupilMath.GeometryCVectorial(20, 0.1, 0.6, 1.5, 1.4)

    pf = geo.createFromZernike(1.0,  [[1.3, -1, 3], [1.3, -2, 2]])    
    z_vals = numpy.linspace(-1.0,1.0,10)

    gsf = geo.gaussianScalingFactor(1.8)
    psf_py = geo.pfToPSF(pf, z_vals, scaling_factor = gsf)
    psf_c = geo_c.pfToPSF(pf, z_vals, scaling_factor = gsf)
    
    
if (__name__ == "__main__"):
    test_pupil_math_1()
    test_pupil_math_2()
    test_pupil_math_3()

    
