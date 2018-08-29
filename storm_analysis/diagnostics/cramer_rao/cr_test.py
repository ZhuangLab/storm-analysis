#!/usr/bin/env python
"""
Test(s) of Cramer-Rao bounds calculations.

Hazen 10/17
"""
import math
import numpy

import storm_analysis.psf_fft.cramer_rao as psfFFTCramerRao
import storm_analysis.pupilfn.cramer_rao as pupilFnCramerRao
import storm_analysis.spliner.cramer_rao as splinerCramerRao

import storm_analysis.diagnostics.cramer_rao.settings as settings


def crTest():
    cr_psf_fft = psfFFTCramerRao.CRPSFFn(psf_filename = "psf_fft.psf",
                                         pixel_size = settings.pixel_size)

    cr_pupil_fn = pupilFnCramerRao.CRPupilFn(psf_filename = "pupilfn.pfn",
                                             pixel_size = settings.pixel_size,
                                             zmax = settings.spline_z_range,
                                             zmin = -settings.spline_z_range)

    cr_spline = splinerCramerRao.CRSplineToPSF3D(psf_filename = "psf.spline",
                                                 pixel_size = settings.pixel_size)

    cr_objects = [cr_spline, cr_psf_fft, cr_pupil_fn]
    
    z_vals = numpy.arange(-settings.test_z_range,
                          settings.test_z_range + 0.5 *settings.test_z_step,
                          settings.test_z_step)

    for i in range(z_vals.size):
        z = z_vals[i]
        print("{0:.3f}".format(z))
        for cro in cr_objects:
            crbs = map(math.sqrt, splinerCramerRao.calcCRBound3D(cro, 20, 1000, z))
            print("{0:.1f} {1:.2f} {2:.2f} {3:.2f} {4:.4f}".format(*crbs))
        print()


if (__name__ == "__main__"):
    crTest()
