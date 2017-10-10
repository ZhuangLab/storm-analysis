#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.simulator.pupil_math as pupilMath
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

    psf_c = pf_c.getPSF()
    psf_py = pupilMath.intensity(pupilMath.toRealSpace(pf))

    assert (numpy.max(numpy.abs(psf_c - psf_py))) < 1.0e-10

    if False:
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_1.tif")) as tf:
            tf.save(psf_c.astype(numpy.float32))
            tf.save(psf_py.astype(numpy.float32))

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

    pf_c.translate(dx, dy, -dz)
    psf_c = pf_c.getPSF()

    defocused = geo.changeFocus(pf, dz)
    translated = geo.translatePf(defocused, dx, dy)
    psf_py = pupilMath.intensity(pupilMath.toRealSpace(translated))

    assert (numpy.max(numpy.abs(psf_c - psf_py))) < 1.0e-10

    if False:
        with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_2.tif")) as tf:
            tf.save(psf_c.astype(numpy.float32))
            tf.save(psf_py.astype(numpy.float32))

    pf_c.cleanup()


if (__name__ == "__main__"):
    test_pupilfn_1()
    test_pupilfn_2()

