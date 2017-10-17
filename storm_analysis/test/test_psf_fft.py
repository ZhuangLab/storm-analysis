#!/usr/bin/env python
import numpy
import tifffile

import storm_analysis
import storm_analysis.psf_fft.psf_fft_c as psfFFTC
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
    pfft.translate(-dy, -dx, dz*(pf_psf.shape[0] - 1)/0.8)
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
    
    if True:
        print(numpy.max(numpy.abs(dx_exact - dx_calc)))
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_psf_fft3.tif")) as tf:
            tf.save(dx_exact.astype(numpy.float32))
            tf.save(dx_calc.astype(numpy.float32))

    assert (numpy.max(numpy.abs(dx_exact - dx_calc))) < 1.0e-6

    pfft.cleanup()    
    
            
if (__name__ == "__main__"):
#    test_psf_fft1()
#    test_psf_fft2()
    test_psf_fft3()
