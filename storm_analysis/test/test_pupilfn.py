#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.simulator.pupil_math as pupilMath
import storm_analysis.pupilfn.pupil_function_c as pfFnC

import tifffile

def test_pupilfn_1():
    geo = pupilMath.Geometry(20, 0.1, 0.6, 1.5, 1.4)
    pf = geo.createFromZernike(1.0, [[1.3, 2, 2]])

    pf_c = pfFnC.PupilFunction(geometry = geo)
    pf_c.setPF(pf)

    psf_c = pf_c.getPSF()
    psf_py = pupilMath.intensity(pupilMath.toRealSpace(pf))

    with tifffile.TiffWriter(storm_analysis.getPathOutputTest("test_pupilfn_1.tif")) as tf:
        tf.save(psf_c.astype(numpy.float32))
        tf.save(psf_py.astype(numpy.float32))

    pf_c.cleanup()



if (__name__ == "__main__"):
    test_pupilfn_1()
