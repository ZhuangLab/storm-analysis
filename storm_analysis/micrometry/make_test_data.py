#!/usr/bin/env python
"""
Generate test data sets for micrometry, two HDF5 localization
files that share some of their localizations.

The shared localizations in the second file are those of the first
file rotated and scaled about the center of the field, optionally
mirrored in x, then shifted. Localizations that this moves outside
the field are dropped.

Hazen 07/17
"""
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py


def makeTestData(locs1, locs2, total, match, im_size = 512, rotation = 0.0, scale = 1.0, dx = 0.0, dy = 0.0,
                 mirror = False, localization_error = 0.0, seed = None):
    """
    locs1 - The name of the first HDF5 file.
    locs2 - The name of the second HDF5 file.
    total - The number of localizations in each file (before dropping
            those outside the field).
    match - The fraction of localizations that are in both files.
    im_size - The size of the (square) field in pixels.
    rotation - The rotation of the second field in degrees.
    scale - The scale of the second field.
    dx, dy - The shift of the second field in pixels.
    mirror - Mirror the second field in x.
    localization_error - The localization error in pixels, x and y.
    seed - Random number generator seed.
    """
    rng = numpy.random.default_rng(seed)

    # Create matching points.
    total_match = int(total * match)
    total_noise = total - total_match
    m_xy = rng.uniform(high = im_size, size = (total_match, 2))

    # Transform the matching points into the second field.
    angle = numpy.deg2rad(rotation)
    matrix = scale * numpy.array([[numpy.cos(angle), -numpy.sin(angle)],
                                  [numpy.sin(angle), numpy.cos(angle)]])
    if mirror:
        matrix = numpy.matmul(matrix, numpy.diag([-1.0, 1.0]))
    center = numpy.array([0.5 * im_size, 0.5 * im_size])
    m2_xy = numpy.matmul(m_xy - center, matrix.T) + center + numpy.array([dx, dy])

    # Add noise points and localization error.
    xy1 = numpy.concatenate((m_xy, rng.uniform(high = im_size, size = (total_noise, 2))))
    xy2 = numpy.concatenate((m2_xy, rng.uniform(high = im_size, size = (total_noise, 2))))
    for [filename, xy] in [[locs1, xy1], [locs2, xy2]]:
        xy = xy + rng.normal(scale = localization_error, size = xy.shape)
        xy = xy[numpy.all((xy >= 0.0) & (xy < im_size), axis = 1)]
        xy = rng.permutation(xy)

        with saH5Py.SAH5Py(filename, is_existing = False, overwrite = True) as h5:
            h5.setMovieInformation(im_size, im_size, 1, "")
            h5.addLocalizations({"x" : xy[:,0], "y" : xy[:,1]}, 0)


if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = 'Generate test data for Micrometry')

    parser.add_argument('--locs1', dest='locs1', type=str, required=False, default="locs1.hdf5",
                        help = "The name of the first localizations file, default is locs1.hdf5.")
    parser.add_argument('--locs2', dest='locs2', type=str, required=False, default="locs2.hdf5",
                        help = "The name of the second localizations file, default is locs2.hdf5.")
    parser.add_argument('--total', dest='total', type=int, required=True,
                        help = "The total number of localizations.")
    parser.add_argument('--match', dest='match', type=float, required=True,
                        help = "The fraction of localizations that match.")
    parser.add_argument('--im_size', dest='im_size', type=int, required=False, default=512,
                        help = "The size of the field in pixels, default is 512.")
    parser.add_argument('--rotation', dest='rotation', type=float, required=False, default=0.0,
                        help = "Rotation of the second field in degrees, default is 0.0.")
    parser.add_argument('--scale', dest='scale', type=float, required=False, default=1.0,
                        help = "Scale of the second field, default is 1.0.")
    parser.add_argument('--dx', dest='dx', type=float, required=False, default=0.0,
                        help = "X shift of the second field in pixels, default is 0.0.")
    parser.add_argument('--dy', dest='dy', type=float, required=False, default=0.0,
                        help = "Y shift of the second field in pixels, default is 0.0.")
    parser.add_argument('--mirror', dest='mirror', action='store_true', default=False,
                        help = "Mirror the second field in x.")
    parser.add_argument('--localization_error', dest='localization_error', type=float, required=False, default=0.0,
                        help = "Localization error in pixels, default is 0.0.")
    parser.add_argument('--seed', dest='seed', type=int, required=False, default=None,
                        help = "Random number generator seed.")

    args = parser.parse_args()

    makeTestData(args.locs1,
                 args.locs2,
                 args.total,
                 args.match,
                 im_size = args.im_size,
                 rotation = args.rotation,
                 scale = args.scale,
                 dx = args.dx,
                 dy = args.dy,
                 mirror = args.mirror,
                 localization_error = args.localization_error,
                 seed = args.seed)
