#!/usr/bin/env python
"""
Verify that the shape of the Cramer-Rao object PSF is correct.

Hazen 10/17
"""
import numpy
import tifffile

import storm_analysis.psf_fft.cramer_rao as psfFFTCramerRao
import storm_analysis.pupilfn.cramer_rao as pupilFnCramerRao
import storm_analysis.spliner.cramer_rao as splinerCramerRao

import storm_analysis.diagnostics.cramer_rao.settings as settings


def ndiff(psf1, psf2):
    return numpy.amax(numpy.abs(psf1 - psf2))/numpy.amax(psf2)

def psfTest():    
    cr_psf_fft = psfFFTCramerRao.CRPSFFn(psf_filename = "psf_fft.psf",
                                         pixel_size = settings.pixel_size)

    cr_pupil_fn = pupilFnCramerRao.CRPupilFn(psf_filename = "pupilfn.pfn",
                                             pixel_size = settings.pixel_size,
                                             zmax = settings.spline_z_range,
                                             zmin = -settings.spline_z_range)

    cr_spline = splinerCramerRao.CRSplineToPSF3D(psf_filename = "psf.spline",
                                                 pixel_size = settings.pixel_size)

    z_vals = numpy.arange(-settings.test_z_range,
                          settings.test_z_range + 0.5 *settings.test_z_step,
                          settings.test_z_step) * 1.0e+3

    if True:
        with tifffile.TiffWriter("psfs.tif") as tf:
            for i in range(z_vals.size):
                z = z_vals[i]
            
                # These two are the same size.
                psf_fft_dx = cr_psf_fft.getDx(z)/cr_psf_fft.delta_xy
                psf_fft_dy = cr_psf_fft.getDy(z)/cr_psf_fft.delta_xy
                psf_fft_dz = cr_psf_fft.getDz(z)/cr_psf_fft.delta_z
                psf_fft = cr_psf_fft.getPSF(z)

                psf_pfn_dx = cr_pupil_fn.getDx(z)/cr_pupil_fn.delta_xy
                psf_pfn_dy = cr_pupil_fn.getDy(z)/cr_pupil_fn.delta_xy
                psf_pfn_dz = cr_pupil_fn.getDz(z)/cr_pupil_fn.delta_z
                psf_pfn = cr_pupil_fn.getPSF(z)

                psf_sp_dx = cr_spline.getDx(z)/cr_spline.delta_xy
                psf_sp_dy = cr_spline.getDy(z)/cr_spline.delta_xy
                psf_sp_dz = cr_spline.getDz(z)/cr_spline.delta_z
                psf_sp = cr_spline.getPSF(z)

                # Check differences:
                print("{0:.3f}".format(z))
                print("dx {0:.3f} {1:.3f}".format(ndiff(psf_fft_dx, psf_sp_dx), ndiff(psf_pfn_dx, psf_sp_dx)))
                print("dy {0:.3f} {1:.3f}".format(ndiff(psf_fft_dy, psf_sp_dy), ndiff(psf_pfn_dy, psf_sp_dy)))
                print("dz {0:.3f} {1:.3f}".format(ndiff(psf_fft_dz, psf_sp_dz), ndiff(psf_pfn_dz, psf_sp_dz)))
                print("psfs {0:.3f} {1:.3f}".format(ndiff(psf_fft, psf_sp), ndiff(psf_pfn, psf_sp)))
                print("")

                # Make a composite image of the 3 different calculations.
                c_x_size = psf_fft.shape[0]
                c_y_size = psf_fft.shape[1]
                composite = numpy.zeros((c_x_size,c_y_size*3))
                for i, elt in enumerate([psf_sp, psf_fft, psf_pfn]):
                    composite[:,i*c_y_size:(i+1)*c_y_size] = elt
            
                tf.save(composite.astype(numpy.float32))

                
if (__name__ == "__main__"):
    psfTest()
