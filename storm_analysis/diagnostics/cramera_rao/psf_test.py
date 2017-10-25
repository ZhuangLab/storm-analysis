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

import settings


def ndiff(psf1, psf2):
    return numpy.amax(numpy.abs(psf1 - psf2))/numpy.amax(psf2)

    
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
                      settings.test_z_step)

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

            # This spline is 1 pixel smaller.
            psf_small = cr_spline.getDx(z)/cr_spline.delta_xy
            psf_sp_dx = numpy.zeros((psf_small.shape[0]+1, psf_small.shape[1]+1))
            psf_sp_dx[1:,1:] = psf_small

            psf_small = cr_spline.getDy(z)/cr_spline.delta_xy
            psf_sp_dy = numpy.zeros((psf_small.shape[0]+1, psf_small.shape[1]+1))
            psf_sp_dy[1:,1:] = psf_small

            psf_small = cr_spline.getDz(z)/cr_spline.delta_z
            psf_sp_dz = numpy.zeros((psf_small.shape[0]+1, psf_small.shape[1]+1))
            psf_sp_dz[1:,1:] = psf_small

            psf_small = cr_spline.getPSF(z)
            psf_sp = numpy.zeros((psf_small.shape[0]+1, psf_small.shape[1]+1))
            psf_sp[1:,1:] = psf_small

            # Check differences:
            print(z)
            print("dx")
            print(ndiff(psf_fft_dx, psf_sp_dx), ndiff(psf_pfn_dx, psf_sp_dx))
            print("dy")
            print(ndiff(psf_fft_dy, psf_sp_dy), ndiff(psf_pfn_dy, psf_sp_dy))
            print("dz")
            print(ndiff(psf_fft_dz, psf_sp_dz), ndiff(psf_pfn_dz, psf_sp_dz))
            print("pfs")
            print(ndiff(psf_fft, psf_sp), ndiff(psf_pfn, psf_sp))
            print("")
            
            tf.save(psf_fft_dz.astype(numpy.float32))
            tf.save(psf_pfn_dz.astype(numpy.float32))
            tf.save(psf_sp_dz.astype(numpy.float32))

