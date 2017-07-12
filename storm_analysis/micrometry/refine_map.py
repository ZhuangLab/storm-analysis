#!/usr/bin/env python
"""
Given two lists of localizations and the 'first guess' from
micrometry this creates a refined mapping.

FIXME: Add degree > 1 mapping capabilities.

Hazen 07/17
"""

import numpy
import pickle
import scipy
import scipy.spatial

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3

import storm_analysis.micrometry.micrometry as micrometry


def calculateTransform(x1, y1, x2, y2):
    """
    Given control points, return the best transform.
    """
    assert(x1.size == x2.size)
    
    m = numpy.ones((x2.size, 3))
    m[:,1] = x2
    m[:,2] = y2

    return [numpy.linalg.lstsq(m, x1)[0],
            numpy.linalg.lstsq(m, y1)[0]]

    
def findControlPoints(kd1, kd2, transform, max_distance = 1.0):
    """
    Returns matching points (within max_distance) in 
    the list [x1, y1, x2, y2]
    """
    # Transform to reference frame.
    [x2, y2] = micrometry.applyTransform(kd2, transform)
    p2 = numpy.stack((x2, y2), axis = -1)

    # Find distance to nearest point in the reference.
    [dist, index] = kd1.query(p2, distance_upper_bound = max_distance)
    mask = (dist != numpy.inf)

    # Extract 'reference' control points.
    index = index[mask]
    x1 = numpy.zeros(numpy.count_nonzero(mask))
    y1 = numpy.zeros(numpy.count_nonzero(mask))
    x1 = kd1.data[index,0].copy()
    y1 = kd1.data[index,1].copy()
    
    # Extract 'other' control points.
    x2 = kd2.data[mask,0].copy()
    y2 = kd2.data[mask,1].copy()

    print(x1.size, "control points found.")
    
    return [x1, y1, x2, y2]


def refineTransform(kd1, kd2, transform):
    """
    Refines the transform.

    kd1 and kd2 are scipy.spatial.KDTree objects.
    transform is an initial guess for the transform.
    
    Returns the refined transform and it's inverse.
    """
    # Update transform.
    transform = calculateTransform(*findControlPoints(kd1, kd2, transform))

    # Find control points.
    [x1, y1, x2, y2] = findControlPoints(kd1, kd2, transform)

    # '0' is the reference, '1' is the other.
    tr_1_to_0 = calculateTransform(x1, y1, x2, y2)
    tr_0_to_1 = calculateTransform(x2, y2, x1, y1)

    return [tr_1_to_0, tr_0_to_1]
    

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

    args = parser.parse_args()

    # Load locs1.
    i3_data1 = readinsight3.loadI3File(args.locs1)
    kd1 = scipy.spatial.KDTree(numpy.stack((i3_data1['xc'], i3_data1['yc']), axis = -1))

    # Load locs2.
    i3_data2 = readinsight3.loadI3File(args.locs2)
    kd2 = scipy.spatial.KDTree(numpy.stack((i3_data2['xc'], i3_data2['yc']), axis = -1))

    # Load 'first guess' transform.
    with open(args.mm_map, 'rb') as fp:
        mp_transform = pickle.load(fp)

    # Refine.
    [tr_1_to_0, tr_0_to_1] = refineTransform(kd1, kd2, [mp_transform["1_0_x"], mp_transform["1_0_y"]])

    # Save results.
    mapping = {"1_0_x" : tr_1_to_0[0],
               "1_0_y" : tr_1_to_0[1],
               "0_1_x" : tr_0_to_1[0],
               "0_1_y" : tr_0_to_1[1]}

    with open(args.results, 'wb') as fp:
        pickle.dump(mapping, fp)
