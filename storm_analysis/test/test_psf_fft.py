#!/usr/bin/env python
import numpy
import tifffile

import storm_analysis
import storm_analysis.psf_fft.psf_fft_c as psfFFTC
import storm_analysis.psf_fft.psf_fft_py as psfFFTPy
import storm_analysis.simulator.pupil_math as pupilMath


def makePSFAndPF(zmin, zmax, zstep):
    """
    Creates the PSF and PF used for testing.
    """
    size = 20
    geo = pupilMath.Geometry(size, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0, [[1.3, 2, 2]])

    z_values = numpy.arange(zmin, zmax + 0.5*zstep, zstep)
    
    psf = numpy.zeros((z_values.size, size, size))
    for i, z in enumerate(z_values):
        defocused = geo.changeFocus(pf, z)
        psf[i,:,:] = pupilMath.intensity(pupilMath.toRealSpace(defocused))

    return [psf, geo, pf]
    
    
def test_psf_fft1():
    """
    Test untranslated PSF calculation.
    """
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)

    pfft = psfFFTC.PSFFFT(pf_psf)
    psf_fft = pfft.getPSF()
    psf_pf = pupilMath.intensity(pupilMath.toRealSpace(pf))
    
    if False:
        print(numpy.max(numpy.abs(psf_fft - psf_pf)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft1.tif")) as tf:
            tf.save(psf_fft.astype(numpy.float32))
            tf.save(psf_pf.astype(numpy.float32))

    assert (numpy.max(numpy.abs(psf_fft - psf_pf))) < 1.0e-10

    pfft.cleanup()

def test_psf_fft2():
    """
    Test translated PSF calculation.
    """
    dx = 0.5
    dy = 0.25
    dz = 0.2
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)
    
    pfft = psfFFTC.PSFFFT(pf_psf)
    pfft.translate(dy, dx, dz*(pf_psf.shape[0] - 1)/0.8)
    psf_fft = pfft.getPSF()
    
    defocused = geo.changeFocus(pf, dz)
    translated = geo.translatePf(defocused, dx, dy)
    psf_pf = pupilMath.intensity(pupilMath.toRealSpace(translated))
    
    if False:
        print(numpy.max(numpy.abs(psf_fft - psf_pf)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft2.tif")) as tf:
            tf.save(psf_fft.astype(numpy.float32))
            tf.save(psf_pf.astype(numpy.float32))

    assert (numpy.max(numpy.abs(psf_fft - psf_pf))) < 1.0e-10

    pfft.cleanup()

def test_psf_fft3():
    """
    Test PSF dx calculation.
    """
    dx = 0.2
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)
    
    pfft = psfFFTC.PSFFFT(pf_psf)
    pfft.translate(dx, 0.0, 0.0)
    dx_exact = pfft.getPSFdx()

    psf = pfft.getPSF()
    pfft.translate(dx + 1.0e-6, 0.0, 0.0)
    dx_calc = (pfft.getPSF() - psf)/1.0e-6
    
    if False:
        print(numpy.max(numpy.abs(dx_exact - dx_calc)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft3.tif")) as tf:
            tf.save(dx_exact.astype(numpy.float32))
            tf.save(dx_calc.astype(numpy.float32))

    assert (numpy.max(numpy.abs(dx_exact - dx_calc))) < 1.0e-6

    pfft.cleanup()

def test_psf_fft4():
    """
    Test PSF dy calculation.
    """
    dy = 0.0
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)
    
    pfft = psfFFTC.PSFFFT(pf_psf)
    pfft.translate(0.0, dy, 0.0)
    dy_exact = pfft.getPSFdy()

    psf = pfft.getPSF()
    pfft.translate(0.0, dy + 1.0e-6, 0.0)
    dy_calc = (pfft.getPSF() - psf)/1.0e-6
    
    if False:
        print(numpy.max(numpy.abs(dy_exact - dy_calc)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft4.tif")) as tf:
            tf.save(dy_exact.astype(numpy.float32))
            tf.save(dy_calc.astype(numpy.float32))

    assert (numpy.max(numpy.abs(dy_exact - dy_calc))) < 1.0e-6

    pfft.cleanup()

def test_psf_fft5():
    """
    Test PSF dz calculation.
    """
    dz = 0.0
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)
    
    pfft = psfFFTC.PSFFFT(pf_psf)
    pfft.translate(0.0, 0.0, dz)
    dz_exact = pfft.getPSFdz()

    psf = pfft.getPSF()
    pfft.translate(0.0, 0.0, dz + 1.0e-6)
    dz_calc = (pfft.getPSF() - psf)/1.0e-6
    
    if False:
        print(numpy.max(numpy.abs(dz_exact - dz_calc)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft5.tif")) as tf:
            tf.save(dz_exact.astype(numpy.float32))
            tf.save(dz_calc.astype(numpy.float32))

    assert (numpy.max(numpy.abs(dz_exact - dz_calc))) < 1.0e-6

    pfft.cleanup()    

def test_psf_fft6():
    """
    Test against the Python version, no translation.
    """
    dz = 0.0
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)
    
    pfft_c = psfFFTC.PSFFFT(pf_psf)
    pfft_py = psfFFTPy.PSFFFT(pf_psf)

    psf_c =pfft_c.getPSF()
    psf_py = pfft_py.getPSF()

    if False:
        print(numpy.max(numpy.abs(psf_c - psf_py)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft6.tif")) as tf:
            tf.save(psf_c.astype(numpy.float32))
            tf.save(psf_py.astype(numpy.float32))

    assert (numpy.max(numpy.abs(psf_c - psf_py))) < 1.0e-6

    pfft_c.cleanup()

def test_psf_fft7():
    """
    Test against the Python version, translation.
    """
    dx = 0.5
    dy = 0.25
    dz = 0.2
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)
    
    pfft_c = psfFFTC.PSFFFT(pf_psf)
    pfft_py = psfFFTPy.PSFFFT(pf_psf)

    pfft_c.translate(dx, dy, dz)
    pfft_py.translate(dx, dy, dz)

    psf_c =pfft_c.getPSF()
    psf_py = pfft_py.getPSF()

    if False:
        print(numpy.max(numpy.abs(psf_c - psf_py)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft7.tif")) as tf:
            tf.save(psf_c.astype(numpy.float32))
            tf.save(psf_py.astype(numpy.float32))

    assert (numpy.max(numpy.abs(psf_c - psf_py))) < 1.0e-6

    pfft_c.cleanup()
    
def test_psf_fft8():
    """
    Test against the Python version, dx.
    """
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)
    
    pfft_c = psfFFTC.PSFFFT(pf_psf)
    pfft_py = psfFFTPy.PSFFFT(pf_psf)

    psf_dx_c = pfft_c.getPSFdx()
    psf_dx_py = pfft_py.getPSFdx()

    if False:
        print(numpy.max(numpy.abs(psf_dx_c - psf_dx_py)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft8.tif")) as tf:
            tf.save(psf_dx_c.astype(numpy.float32))
            tf.save(psf_dx_py.astype(numpy.float32))

    assert (numpy.max(numpy.abs(psf_dx_c - psf_dx_py))) < 1.0e-6

    pfft_c.cleanup()

def test_psf_fft9():
    """
    Test against the Python version, dy.
    """
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)
    
    pfft_c = psfFFTC.PSFFFT(pf_psf)
    pfft_py = psfFFTPy.PSFFFT(pf_psf)

    psf_dy_c = pfft_c.getPSFdy()
    psf_dy_py = pfft_py.getPSFdy()

    if False:
        print(numpy.max(numpy.abs(psf_dy_c - psf_dy_py)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft9.tif")) as tf:
            tf.save(psf_dy_c.astype(numpy.float32))
            tf.save(psf_dy_py.astype(numpy.float32))

    assert (numpy.max(numpy.abs(psf_dy_c - psf_dy_py))) < 1.0e-6

    pfft_c.cleanup()

def test_psf_fft10():
    """
    Test against the Python version, dz.
    """
    [pf_psf, geo, pf] = makePSFAndPF(-0.4, 0.4, 0.05)
    
    pfft_c = psfFFTC.PSFFFT(pf_psf)
    pfft_py = psfFFTPy.PSFFFT(pf_psf)

    psf_dz_c = pfft_c.getPSFdz()
    psf_dz_py = pfft_py.getPSFdz()

    if False:
        print(numpy.max(numpy.abs(psf_dz_c - psf_dz_py)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft10.tif")) as tf:
            tf.save(psf_dz_c.astype(numpy.float32))
            tf.save(psf_dz_py.astype(numpy.float32))

    assert (numpy.max(numpy.abs(psf_dz_c - psf_dz_py))) < 1.0e-6

    pfft_c.cleanup()
    
if (__name__ == "__main__"):
    test_psf_fft1()
    test_psf_fft2()
    test_psf_fft3()
    test_psf_fft4()
    test_psf_fft5()
    test_psf_fft6()
    test_psf_fft7()
    test_psf_fft8()
    test_psf_fft9()
    test_psf_fft10()
    
