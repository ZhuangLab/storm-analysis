#!/usr/bin/env python
"""
Creates a localization file with a random distribution of emitters.

Hazen 09/17
"""
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py


def emittersUniformRandom(h5_name, density, margin, sx, sy, zrange, seed = 0, sigma = 1.5):
    """
    h5_name - The name of the HDF5 file to save the emitter locations, etc.
    density - Localizations per pixel squared.
    margin - The margin in pixels around the edge of the image.
    sx - The image size in Y in pixels.
    sy - The image size in Y in pixels.
    zrange - Range for z values in microns, -zrange to zrange.
    seed - A seed for the random number generator, default is 0.
    sigma - The sigma for the localizatoins, default is 1.5 pixels.
    """
    
    # Set random number generator seed.
    if seed is not None:
        numpy.random.seed(seed)

    # Calculate the size of area covered by the localizations.
    size_x = sx - 2*margin
    size_y = sy - 2*margin

    # Calculate number of localizations.
    n_locs = int(round(size_x*size_y*density))

    # Create localizations.
    peaks = {}
    peaks["x"] = margin + size_x * numpy.random.uniform(size = n_locs)
    peaks["y"] = margin + size_y * numpy.random.uniform(size = n_locs)
    peaks["z"] = -zrange + 2.0 * zrange * numpy.random.uniform(size = n_locs)
    peaks["xsigma"] = sigma*numpy.ones(peaks["x"].size)
    peaks["ysigma"] = sigma*numpy.ones(peaks["y"].size)

    # Save localizations.
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(sx, sy, 1, "")
        h5.addLocalizations(peaks, 0)


if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = "Create a uniform random distribution of emitters for simulations.")
    
    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the HDF5 file to save the emitter locations, etc.")
    parser.add_argument('--density', dest='density', type=float, required=True,
                        help = "Localizations per pixel squared.")
    parser.add_argument('--margin', dest='margin', type=int, required=False, default = 10,
                        help = "The margin in pixels around the edge of the image, default is 10.")
    parser.add_argument('--sx', dest='sx', type=int, required=False, default = 256,
                        help = "The image size in Y in pixels, default is 256.")
    parser.add_argument('--sy', dest='sy', type=int, required=False, default = 256,
                        help = "The image size in Y in pixels, default is 256.")
    parser.add_argument('--zrange', dest='zrange', type=float, required=False, default = 0.0,
                        help = "Range for z values in microns, -zrange to zrange, default is 0.0.")

    args = parser.parse_args()

    emittersUniformRandom(args.hdf5, args.density, args.margin, args.sx, args.sy, args.zrange)



