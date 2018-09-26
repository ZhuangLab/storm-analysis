#!/usr/bin/env python
"""
Create a pupil function file that can be used by pupilfn_analysis for
SMLM analysis. The pupil function is normalized to give a PSF with
a height of maximum height of 1.0.

This is primarily designed for testing the analysis, so it will default
to using the pupil_math.GeometrySim class. See note in simulator.pupil_math
for a more detailed explanation.

Hazen 10/17
"""
import ast
import pickle
import math
import numpy

import storm_analysis.simulator.psf as simPSF
import storm_analysis.simulator.pupil_math as pupilMath


def makePupilFunction(filename, size, pixel_size, zmn, z_offset = 0.0, geo_sim_pf = True):
    """
    filename - The name of the file to save the pupil function in.
    size - The size of the pupil function in pixels.
    pixel_size - pixel size in microns.
    zmn - Zernike coefficients.
    z_offset - Amount to change the focus by in microns.
    geo_sim_pf - Use the 'simulation' PF with 1/2 the pixel size.
    """    
    # This is a requirement of the C library.
    assert ((size%2)==0)
    
    # Physical constants. Note that these match the default values for
    # simulator.psf.PupilFunction().
    #
    wavelength = 1.0e-3 * simPSF.pf_wavelength  # Fluorescence wavelength in microns.
    imm_index = simPSF.pf_refractive_index      # Immersion media index (oil objective).
    NA = simPSF.pf_numerical_aperture           # Numerical aperture of the objective.

    # Create geometry object.
    if geo_sim_pf:
        geo = pupilMath.GeometrySim(size,
                                    pixel_size,
                                    wavelength,
                                    imm_index,
                                    NA)
        
    else:
        geo = pupilMath.Geometry(size,
                                 pixel_size,
                                 wavelength,
                                 imm_index,
                                 NA)

    # Create PF.
    pf = geo.createFromZernike(1.0, zmn)

    # Normalize to have height 1.0.
    psf = pupilMath.intensity(pupilMath.toRealSpace(pf))
    pf = pf * 1.0/math.sqrt(numpy.max(psf))

    # Verify normalization.
    print("Height:", numpy.max(pupilMath.intensity(pupilMath.toRealSpace(pf))))

    # Heh, if zmn is an empty list the pupil function will be perfectly
    # symmetric at z = 0 and the solver will fail because dz = 0. So we
    # solve this we adding a little noise.
    if (len(zmn) == 0):
        print("Plane wave PF detected! Adding noise to break z = 0 symmetry!")
        n_mag = numpy.real(pf) * 1.0e-3
        pf = pf + n_mag * (numpy.random.uniform(size = pf.shape) - 0.5)
    
    # Change focus by z_offset.
    #
    # The convention is that if z_offset + localization z is the final
    # z position, so if localization z is = -z_offset then you will get
    # the PSF at z = 0.
    #
    pf = geo.changeFocus(pf, -z_offset)

    # Pickle and save.
    pfn_dict = {"pf" : pf,
                "pixel_size" : pixel_size,
                "geo_sim_pf" : geo_sim_pf,
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
    parser.add_argument('--sim-pf', dest='geo_sim_pf', type=int, required=False, default = 1,
                        help = "Use PF modified for simulations (0|1), default is 1.")
    parser.add_argument('--zmn', dest='zmn', type=str, required=False, default = "[[1.3, 2, 2]]",
                        help = "The Zernike polynomial coefficient.")
    parser.add_argument('--z-offset', dest='z_offset', type=float, required=False, default = 0.0,
                        help = "Focal plane offset in microns.")

    args = parser.parse_args()

    makePupilFunction(args.filename,
                      args.size,
                      args.pixel_size * 1.0e-3,
                      ast.literal_eval(args.zmn),
                      z_offset = args.z_offset,
                      geo_sim_pf = (args.geo_sim_pf != 0))
