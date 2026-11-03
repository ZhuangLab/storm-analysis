#!/usr/bin/env python
"""
Given two lists of localizations and a mapping between them, this
creates a refined mapping.

micrometry refines the mapping that it finds, so this is for when
the mapping is already known approximately, for example to update
it with beads imaged on another day.

FIXME: Add degree > 1 mapping capabilities.

Hazen 07/17
"""
import pickle
import scipy
import scipy.spatial

import storm_analysis.micrometry.micrometry as micrometry


def refineMap(locs1, locs2, mm_map, results, localization_error = 0.5):
    """
    locs1 - The name of the 'reference' HDF5 format localizations file.
    locs2 - The name of the 'other' HDF5 format localizations file.
    mm_map - The name of the mapping file to refine.
    results - The name of the file to save the refined mapping in.
    localization_error - The localization error in pixels.

    Returns the number of localizations that were paired.
    """
    kd_ref = scipy.spatial.KDTree(micrometry.loadLocalizations(locs1))
    other_points = micrometry.loadLocalizations(locs2)

    with open(mm_map, 'rb') as fp:
        mapping = pickle.load(fp)

    [transform, pairs] = micrometry.refineTransform(kd_ref,
                                                    other_points,
                                                    [mapping["1_0_x"], mapping["1_0_y"]],
                                                    micrometry.pairingRadius(localization_error))

    print(pairs[0].size, "control points found.")
    if (pairs[0].size < 4):
        print("Not enough control points, the mapping was not saved.")
    else:
        reverse = micrometry.reverseTransform(kd_ref, other_points, transform, pairs)
        micrometry.saveMapping(list(transform) + list(reverse), results)

    return pairs[0].size


if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = 'Refine a mapping.')

    parser.add_argument('--locs1', dest='locs1', type=str, required=True,
                        help = "The name of the 'reference' localizations file")
    parser.add_argument('--locs2', dest='locs2', type=str, required=True,
                        help = "The name of the 'other' localizations file")
    parser.add_argument('--mm_map', dest='mm_map', type=str, required=True,
                        help = "The name of the micrometry map file.")
    parser.add_argument('--results', dest='results', type=str, required=True,
                        help = "The name of the file to save the updated mapping in.")
    parser.add_argument('--localization_error', dest='localization_error', type=float, required=False, default=0.5,
                        help = "Localization error (pixels) in x and y, default is 0.5.")

    args = parser.parse_args()

    refineMap(args.locs1,
              args.locs2,
              args.mm_map,
              args.results,
              localization_error = args.localization_error)
