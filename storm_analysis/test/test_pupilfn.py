#!/usr/bin/env python
import math
import numpy
import pickle

import storm_analysis

import storm_analysis.simulator.pupil_math as pupilMath
import storm_analysis.pupilfn.make_pupil_fn as makePupilFn
import storm_analysis.pupilfn.pupil_function_c as pfFnC

import tifffile

def test_pupilfn_1():
    """
    Test that the C and the Python library agree on the calculation
    of the untranslated PSF.
    """
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0, [[1.3, 2, 2]])

    pf_c = pfFnC.PupilFunction(geometry = geo)
    pf_c.setPF(pf)

    psf_c = pupilMath.intensity(pf_c.getPSF())
    psf_py = pupilMath.intensity(pupilMath.toRealSpace(pf))
    
    if False:
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_1.tif")) as tf:
            tf.save(psf_c.astype(numpy.float32))
            tf.save(psf_py.astype(numpy.float32))

    assert numpy.allclose(psf_c, psf_py)
    
    pf_c.cleanup()

def test_pupilfn_2():
    """
    Test PF translation.
    """
    dx = 0.5
    dy = 0.25
    dz = 0.2
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0, [[1.3, 2, 2]])

    pf_c = pfFnC.PupilFunction(geometry = geo)
    pf_c.setPF(pf)

    pf_c.translate(dx, dy, dz)
    psf_c = pupilMath.intensity(pf_c.getPSF())

    defocused = geo.changeFocus(pf, dz)
    translated = geo.translatePf(defocused, dx, dy)
    psf_py = pupilMath.intensity(pupilMath.toRealSpace(translated))

    if False:
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_2.tif")) as tf:
            tf.save(psf_c.astype(numpy.float32))
            tf.save(psf_py.astype(numpy.float32))

    assert numpy.allclose(psf_c, psf_py)
            
    pf_c.cleanup()

def test_pupilfn_3():
    """
    Test PF X derivative (C library).
    """
    dx = 1.0e-6
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0, [[1.3, 2, 2]])

    pf_c = pfFnC.PupilFunction(geometry = geo)
    pf_c.setPF(pf)
    
    # Calculate derivative of magnitude as a function of x.
    psf_c = pf_c.getPSF()
    psf_c_dx = pf_c.getPSFdx()
    mag_dx_calc = 2.0 * (numpy.real(psf_c)*numpy.real(psf_c_dx) + numpy.imag(psf_c)*numpy.imag(psf_c_dx))

    # Estimate derivative using (f(x+dx) - f(x))/dx
    mag = pupilMath.intensity(psf_c)
    pf_c.translate(dx,0.0,0.0)
    mag_dx_est = (pupilMath.intensity(pf_c.getPSF()) - mag)/dx
                
    if False:
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_3.tif")) as tf:
            #tf.save(mag.astype(numpy.float32))
            tf.save(mag_dx_calc.astype(numpy.float32))
            tf.save(mag_dx_est.astype(numpy.float32))
            tf.save(numpy.abs(mag_dx_calc - mag_dx_est).astype(numpy.float32))

    assert numpy.allclose(mag_dx_calc, mag_dx_est, atol = 1.0e-6)
    
    pf_c.cleanup()

def test_pupilfn_4():
    """
    Test PF X derivative (Python library).
    """
    dx = 1.0e-6
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0, [[1.3, 2, 2]])
    
    # Calculate derivative of magnitude as a function of x.
    psf_py = pupilMath.toRealSpace(pf)
    psf_py_dx = pupilMath.toRealSpace(geo.dx(pf))
    mag_dx_calc = 2.0 * (numpy.real(psf_py)*numpy.real(psf_py_dx) + numpy.imag(psf_py)*numpy.imag(psf_py_dx))

    # Estimate derivative using (f(x+dx) - f(x))/dx
    mag = pupilMath.intensity(psf_py)
    translated = geo.translatePf(pf, dx, 0.0)
    mag_dx_est = (pupilMath.intensity(pupilMath.toRealSpace(translated)) - mag)/dx
        
    if False:
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_4.tif")) as tf:
            #tf.save(mag.astype(numpy.float32))
            tf.save(mag_dx_calc.astype(numpy.float32))
            tf.save(mag_dx_est.astype(numpy.float32))
            tf.save(numpy.abs(mag_dx_calc - mag_dx_est).astype(numpy.float32))

    assert numpy.allclose(mag_dx_calc, mag_dx_est, atol = 1.0e-6)

def test_pupilfn_5():
    """
    Test PF Y derivative (C library).
    """
    dy = 1.0e-6
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0, [[1.3, 2, 2]])

    pf_c = pfFnC.PupilFunction(geometry = geo)
    pf_c.setPF(pf)
    
    # Calculate derivative of magnitude as a function of y.
    psf_c = pf_c.getPSF()
    psf_c_dy = pf_c.getPSFdy()
    mag_dy_calc = 2.0 * (numpy.real(psf_c)*numpy.real(psf_c_dy) + numpy.imag(psf_c)*numpy.imag(psf_c_dy))

    # Estimate derivative using (f(y+dy) - f(y))/dy
    mag = pupilMath.intensity(psf_c)
    pf_c.translate(0.0,dy,0.0)
    mag_dy_est = (pupilMath.intensity(pf_c.getPSF()) - mag)/dy
                
    if False:
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_5.tif")) as tf:
            #tf.save(mag.astype(numpy.float32))
            tf.save(mag_dy_calc.astype(numpy.float32))
            tf.save(mag_dy_est.astype(numpy.float32))
            tf.save(numpy.abs(mag_dy_calc - mag_dy_est).astype(numpy.float32))

    assert numpy.allclose(mag_dy_calc, mag_dy_est, atol = 1.0e-6)
    
    pf_c.cleanup()

def test_pupilfn_6():
    """
    Test PF Z derivative (C library).
    """
    dz = 1.0e-6
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0, [[1.3, 2, 2]])

    pf_c = pfFnC.PupilFunction(geometry = geo)
    pf_c.setPF(pf)
    
    # Calculate derivative of magnitude as a function of z.
    psf_c = pf_c.getPSF()
    psf_c_dz = pf_c.getPSFdz()
    mag_dz_calc = 2.0 * (numpy.real(psf_c)*numpy.real(psf_c_dz) + numpy.imag(psf_c)*numpy.imag(psf_c_dz))

    # Estimate derivative using (f(z+dz) - f(z))/dz
    mag = pupilMath.intensity(psf_c)
    pf_c.translate(0.0,0.0,dz)
    mag_dz_est = (pupilMath.intensity(pf_c.getPSF()) - mag)/dz
                
    if False:
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_6.tif")) as tf:
            #tf.save(mag.astype(numpy.float32))
            tf.save(mag_dz_calc.astype(numpy.float32))
            tf.save(mag_dz_est.astype(numpy.float32))
            tf.save(numpy.abs(mag_dz_calc - mag_dz_est).astype(numpy.float32))

    assert numpy.allclose(mag_dz_calc, mag_dz_est, atol = 1.0e-6)
    
    pf_c.cleanup()    

def test_pupilfn_7():
    """
    Test that PF translation is correct (i.e. independent of size).
    """
    sizes = [10, 20, 40]
    dx = 1.0

    for size in sizes:
        geo = pupilMath.Geometry(size, 0.1, 0.6, 1.5, 1.4)
        pf = geo.createFromZernike(1.0, [[1.3, 2, 2]])

        pf_c = pfFnC.PupilFunction(geometry = geo)
        pf_c.setPF(pf)
        
        psf_untranslated = numpy.roll(pupilMath.intensity(pf_c.getPSF()), 1, axis = 0)
            
        pf_c.translate(dx, 0.0, 0.0)
        psf_translated = pupilMath.intensity(pf_c.getPSF())

        if False:
            with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_7.tif")) as tf:
                tf.save(psf_untranslated.astype(numpy.float32))
                tf.save(psf_translated.astype(numpy.float32))

        assert numpy.allclose(psf_untranslated, psf_translated)
            
        pf_c.cleanup()

def test_pupilfn_8():
    """
    Test that pupilfn.make_pupil_fn.makePupilFunction works as expected.
    """
    pf_size = 30
    zmn = [[1.3, 2, 2]]
    z_offset = -0.3
    
    # Create & save pupil function.
    pf_file = storm_analysis.getPathOutputTest("pf_test.pfn")
    makePupilFn.makePupilFunction(pf_file, pf_size, 0.1, zmn, z_offset = z_offset)

    # Load PF.
    with open(pf_file, "rb") as fp:
        pf_data = pickle.load(fp)
        test_pf = pf_data["pf"]

    # Create comparison PF.
    geo = pupilMath.GeometrySim(pf_size,
                                pf_data["pixel_size"],
                                pf_data["wavelength"],
                                pf_data["immersion_index"],
                                pf_data["numerical_aperture"])
    ref_pf = geo.createFromZernike(1.0, zmn)

    # Normalize reference to also have height 1.0 (at z = 0.0).
    psf = pupilMath.intensity(pupilMath.toRealSpace(ref_pf))
    ref_pf = ref_pf * 1.0/math.sqrt(numpy.max(psf))

    # Test that they are the same.
    for z in [-0.2, -0.1, 0.0, 0.1, 0.2]:
        test_psf = pupilMath.intensity(pupilMath.toRealSpace(geo.changeFocus(test_pf, z)))
        ref_psf = pupilMath.intensity(pupilMath.toRealSpace(geo.changeFocus(ref_pf, z - z_offset)))
        #print(numpy.max(numpy.abs(test_psf - ref_psf)))
        assert numpy.allclose(test_psf, ref_psf)

def test_pupilfn_9():
    """
    Another test PF translation with a less symmetric PSF.
    """
    dx = 0.5
    dy = 0.25
    dz = 0.2
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0,  [[1.3, -1, 3], [1.3, -2, 2]])

    pf_c = pfFnC.PupilFunction(geometry = geo)
    pf_c.setPF(pf)

    pf_c.translate(dx, dy, dz)
    psf_c = pupilMath.intensity(pf_c.getPSF())

    defocused = geo.changeFocus(pf, dz)
    translated = geo.translatePf(defocused, dx, dy)
    psf_py = pupilMath.intensity(pupilMath.toRealSpace(translated))

    if False:
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_2.tif")) as tf:
            tf.save(psf_c.astype(numpy.float32))
            tf.save(psf_py.astype(numpy.float32))

    assert numpy.allclose(psf_c, psf_py)

    pf_c.cleanup()

def test_pupilfn_10():
    """
    Test C library PSF intensity calculation.
    """
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0,  [[1.3, -1, 3], [1.3, -2, 2]])

    pf_c = pfFnC.PupilFunction(geometry = geo)
    pf_c.setPF(pf)

    psf_c = pf_c.getPSFIntensity()
    psf_py = pupilMath.intensity(pupilMath.toRealSpace(pf))

    assert numpy.allclose(psf_c, psf_py)

    pf_c.cleanup()

def test_pupilfn_11():
    """
    Test C library PF Z translation function.
    """
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0,  [[1.3, -1, 3], [1.3, -2, 2]])

    pf_c = pfFnC.PupilFunction(geometry = geo)
    pf_c.setPF(pf)

    for dz in [-0.2, 0.1, 0.0, 0.1, 0.2]:
        pf_c.translateZ(dz)
        psf_c = pf_c.getPSFIntensity()

        defocused = geo.changeFocus(pf, dz)
        psf_py = pupilMath.intensity(pupilMath.toRealSpace(defocused))

        assert numpy.allclose(psf_c, psf_py)

    pf_c.cleanup()
    
            
if (__name__ == "__main__"):
    test_pupilfn_1()
    test_pupilfn_2()
    test_pupilfn_3()
    test_pupilfn_4()
    test_pupilfn_5()
    test_pupilfn_6()
    test_pupilfn_7()
    test_pupilfn_8()
    test_pupilfn_9()
    test_pupilfn_10()
    test_pupilfn_11()
    
    
