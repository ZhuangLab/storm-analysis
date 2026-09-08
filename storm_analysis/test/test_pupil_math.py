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
    
    
def test_pupil_math_aberration_opd():
    """
    aberrationOPD() had two names wrong and could not run at all.

    sin_theta_2 was computed from n1 and n2, neither of which exists, and
    the return statement called self.appylNARestriction(). The first
    raised NameError, and fixing only that moved the failure to an
    AttributeError from the second.

    Snell's law here is between the immersion medium and the sample,
    which is what the sibling aberration() uses:

        sin_theta_2 = (self.imm_index/smp_index)*sin_theta_1

    This geometry deliberately puts most of the grid outside the NA, so
    that the restriction being applied is observable rather than vacuous.
    """
    geo = pupilMath.Geometry(32, 0.05, 0.6, 1.5, 1.2)

    # applyNARestriction() masks on the normalized radius, geo.r > 1.0.
    # geo.r_max is a different quantity, the NA edge in grid units.
    outside = (geo.r > 1.0)
    assert(numpy.count_nonzero(outside) > 0)

    ab = geo.aberrationOPD(1.0, 0.5, 1.33)

    assert(ab.shape == (32, 32))
    assert(numpy.all(numpy.isfinite(ab)))

    # applyNARestriction() zeros everything beyond the NA.
    assert(numpy.all(ab[outside] == 0.0))
    assert(numpy.any(ab[~outside] != 0.0))

    # At zero defocus and a matched sample index there is no aberration,
    # so the function is identically 1 inside the NA.
    flat = geo.aberrationOPD(0.0, 0.0, 1.5)
    assert(numpy.allclose(flat[~outside], 1.0 + 0j))


if (__name__ == "__main__"):
    test_pupil_math_1()
    test_pupil_math_2()
    test_pupil_math_3()
    test_pupil_math_4()
    test_pupil_math_5()
    test_pupil_math_aberration_opd()

    
