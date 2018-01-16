#!/usr/bin/env python
"""
Create a PSF file that can be used for PSF FFT analysis using a 
pupil function.

This is primarily designed for testing the analysis.

Hazen 10/17
"""
import pickle
import math
import numpy
import tifffile

import storm_analysis.simulator.pupil_math as pupilMath


def makePSF(filename, size, pixel_size, zmn, zrange, zstep):
    """
    pixel_size - pixel size in microns.
    zmn - Zernike coefficients.
    zrange - The final zrange will be +- zrange (microns).
    zstep - The z step size in microns.
    """
    #
    # Physical constants. Note that these match the default values for
    # simulator.psf.PupilFunction().
    #
    wavelength = 0.6   # Fluorescence wavelength in microns.
    imm_index = 1.5    # Immersion media index (oil objective).
    NA = 1.4           # Numerical aperture of the objective.

    # Create geometry object.
    geo = pupilMath.Geometry(size,
                             pixel_size,
                             wavelength,
                             imm_index,
                             NA)

    # Create PF.
    pf = geo.createFromZernike(1.0, zmn)

    # Normalize to have height 1.0 (at z = 0.0).
    psf = pupilMath.intensity(pupilMath.toRealSpace(pf))
    pf = pf * 1.0/math.sqrt(numpy.max(psf))
    
    # Verify normalization.
    print("Height:", numpy.max(pupilMath.intensity(pupilMath.toRealSpace(pf))))

    # Create a PSF at each z value.
    z_values = numpy.arange(-zrange, zrange + 0.5*zstep, zstep)
    #print(z_values)
    
    if ((z_values.size%2)==0):
        print("The number of z slices must be an odd number.")
        assert False, "PSF creation failed."
    
    psf = numpy.zeros((z_values.size, size, size))
    for i, z in enumerate(z_values):
        defocused = geo.changeFocus(pf, z)
        psf[i,:,:] = pupilMath.intensity(pupilMath.toRealSpace(defocused))
    
    # Pickle and save.
    psf_dict = {"psf" : psf,
                "pixel_size" : pixel_size,
                "zmax" : 1000.0 * zrange,
                "zmin" : -1000.0 * zrange}

    with open(filename, 'wb') as fp:
        pickle.dump(psf_dict, fp)

    # Also save a .tif version.
    with tifffile.TiffWriter("psf.tif") as tf:
        for i in range(psf.shape[0]):
            tf.save(psf[i,:,:].astype(numpy.float32))
        

if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Make a PSF')

    parser.add_argument('--filename', dest='filename', type=str, required=True,
                        help = "The name of the file to save the pupil function in.")
    parser.add_argument('--size', dest='size', type=int, required=True,
                        help = "The size of the pixel function in pixels.")
    parser.add_argument('--pixel-size', dest='pixel_size', type=float, required=True,
                        help = "The pixel size in nanometers.")
    parser.add_argument('--zrange', dest='zrange', type=float, required=True,
                        help = "The PSF z range in microns.")
    parser.add_argument('--zstep', dest='zstep', type=float, required=False, default = 0.1,
                        help = "The PSF z step size in microns.")

    args = parser.parse_args()

    # FIXME: Make zmn an adjustable parameter.
    makePSF(args.filename,
            args.size,
            args.pixel_size * 1.0e-3,
            [[1.3, 2, 2]],
            args.zrange,
            args.zstep)
