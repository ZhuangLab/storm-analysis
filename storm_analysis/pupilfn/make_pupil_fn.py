#!/usr/bin/env python
"""
Create a pupil function file that can be used by pupilfn_analysis for
SMLM analysis. The pupil function is normalized to give a PSF with
a height of maximum height of 1.0.

This is primarily designed for testing the analysis.

Hazen 10/17
"""
import pickle
import math
import numpy

import storm_analysis.simulator.pupil_math as pupilMath


def makePupilFunction(filename, size, pixel_size, zmn, z_offset = 0.0):
    """
    pixel_size - pixel size in microns.
    zmn - Zernike coefficients.
    z_offset - Amount to change the focus by in microns.
    """

    # This is a requirement of the C library.
    assert ((size%2)==0)
    
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

    # Change focus by z_offset.
    pf = geo.changeFocus(pf, z_offset)

    # Pickle and save.
    pfn_dict = {"pf" : numpy.transpose(pf),
                "pixel_size" : pixel_size,
                "wavelength" : wavelength,
                "immersion_index" : imm_index,
                "numerical_aperture" : NA}

    with open(filename, 'wb') as fp:
        pickle.dump(pfn_dict, fp)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Make a pupil function')

    parser.add_argument('--filename', dest='filename', type=str, required=True,
                        help = "The name of the file to save the pupil function in.")
    parser.add_argument('--size', dest='size', type=int, required=True,
                        help = "The size of the pixel function in pixels.")
    parser.add_argument('--pixel-size', dest='pixel_size', type=float, required=True,
                        help = "The pixel size in nanometers.")
    parser.add_argument('--zmn', dest='zmn', type=str, required=False, default = "[[1.3, 2, 2]]",
                        help = "The Zernike polynomial coefficient.")
    parser.add_argument('--z-offset', dest='z_offset', type=float, required=False, default = 0.0,
                        help = "Focal plane offset in nanometers.")

    args = parser.parse_args()

    makePupilFunction(args.filename,
                      args.size,
                      args.pixel_size * 1.0e-3,
                      eval(args.zmn),
                      z_offset = args.z_offset * 1.0e-3)
