#!/usr/bin/env python
"""
Given two lists of localizations, returns a 'first guess' at the
transform between them. This is degree 1 affine transform. As
such this could fail for images with large differences in field
curvature.

This uses the ideas in this paper, applied to fiducial references
like flourescent beads:

Lang, D., Hogg, D. W., Mierle, K., Blanton, M., & Roweis, S., 2010, 
Astrometry.net: Blind astrometric calibration of arbitrary astronomical 
images, The Astronomical Journal 139, 1782–1800.

Hazen 07/17
"""

import math
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import pickle
import scipy
import scipy.spatial

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3

import storm_analysis.micrometry.quads as quads


def applyTransform(kd, transform):
    tx = transform[0]
    ty = transform[1]
    x = tx[0] + tx[1]*kd.data[:,0] + tx[2]*kd.data[:,1]
    y = ty[0] + ty[1]*kd.data[:,0] + ty[2]*kd.data[:,1]
    return [x, y]


def bgProbability(kd1):
    """
    Returns an estimate of the background, i.e. how likely the
    two data sets aligned by chance. This assumes a uniform
    distribution of points in 'reference' and 'other'.

    FIXME: Should use the image size instead of estimating it.
    """

    # Estimate density of points in 'reference'.
    xmin = numpy.min(kd1.data[:,0])
    xmax = numpy.max(kd1.data[:,0])
    ymin = numpy.min(kd1.data[:,1])
    ymax = numpy.max(kd1.data[:,1])
    density = 1.0/((xmax - xmin)*(ymax - ymin))

    return density


def fgProbability(kd1, kd2, transform):
    """
    Returns an estimate of how likely the transform is correct.
    """
    # Transform 'other' coordinates into the 'reference' frame.
    [x2, y2] = applyTransform(kd2, transform)
    p2 = numpy.stack((x2, y2), axis = -1)

    # Calculate distance to nearest point in 'reference'.
    [dist, index] = kd1.query(p2)

    # Score assuming a localization accuracy of 1 pixel.
    bg_p = bgProbability(kd1)
    fg_p = bg_p + (1.0 - bg_p) * numpy.sum(numpy.exp(-dist*dist*0.5))/float(x2.size)
    return fg_p

    
def makeTreeAndQuads(x, y, min_size = None, max_size = None, max_neighbors = 10):
    """
    Make a KD tree and a list of quads from x, y points.
    """
    kd = scipy.spatial.KDTree(numpy.stack((x, y), axis = -1))
    m_quads = quads.makeQuads(kd,
                              min_size = min_size,
                              max_size = max_size,
                              max_neighbors = max_neighbors)
    return [kd, m_quads]


def makeTreeAndQuadsFromI3File(i3_filename, min_size = None, max_size = None, max_neighbors = 10):
    """
    Make a KD tree and a list of quads from an Insight3 file.

    Note: This file should probably only have localizations for a single frame.
    """
    i3_data = readinsight3.loadI3File(i3_filename)

    # Warning if there is more than 1 frame in the data.
    if (len(numpy.unique(i3_data['fr'])) > 1):
        print("Warning: Localizations in multiple frames detected!")

    return makeTreeAndQuads(i3_data['xc'],
                            i3_data['yc'],
                            min_size = min_size,
                            max_size = max_size,
                            max_neighbors = max_neighbors)


def plotMatch(kd1, kd2, transform, save_as, show = True):
    [x2, y2] = applyTransform(kd2, transform)
    
    fig = pyplot.figure()
    pyplot.scatter(kd1.data[:,0], kd1.data[:,1], facecolors = 'none', edgecolors = 'red', s = 100)
    pyplot.scatter(x2, y2, color = 'green', marker = '+', s = 100)

    legend = pyplot.legend(('reference', 'other'), loc=1)
    pyplot.xlabel("pixels")
    pyplot.ylabel("pixels")

    ax = pyplot.gca()
    ax.set_aspect('equal')

    fig.savefig(save_as)
    
    if show:
        pyplot.show()


if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = 'Micrometry - ...')

    parser.add_argument('--locs1', dest='locs1', type=str, required=True,
                        help = "The name of the 'reference' localizations file")
    parser.add_argument('--locs2', dest='locs2', type=str, required=True,
                        help = "The name of the 'other' localizations file")
    parser.add_argument('--results', dest='results', type=str, required=True,
                        help = "The name of the file to save the transform (if any) in.")    
    parser.add_argument('--min_size', dest='min_size', type=float, required=False, default=5.0,
                        help = "Minimum quad size (pixels), default is 5.0.")
    parser.add_argument('--max_size', dest='max_size', type=float, required=False, default=100.0,
                        help = "Maximum quad size (pixels), default is 100.0.")
    parser.add_argument('--max_neighbors', dest='max_neighbors', type=int, required=False, default=20,
                        help = "Maximum neighbors to search when making quads.")
    parser.add_argument('--tolerance', dest='tolerance', type=float, required=False, default=1.0e-2,
                        help = "Tolerance for matching quads, default is 1.0e-2.")
    parser.add_argument('--no_plots', dest='no_plots', type=bool, required=False, default=False,
                        help = "Don't show plot of the results.")

    args = parser.parse_args()

    print("Making quads for the 'reference' data.")
    [kd1, quads1] = makeTreeAndQuadsFromI3File(args.locs1,
                                               min_size = args.min_size,
                                               max_size = args.max_size,
                                               max_neighbors = args.max_neighbors)
    print("Created", len(quads1), "quads")
    print("")

    print("Making quads for the 'other' data.")
    [kd2, quads2] = makeTreeAndQuadsFromI3File(args.locs2,
                                               min_size = args.min_size,
                                               max_size = args.max_size,
                                               max_neighbors = args.max_neighbors)
    print("Created", len(quads2), "quads")
    print("")
    
    print("Comparing quads.")
    bg_p = bgProbability(kd1)

    #
    # Unlike astrometry.net we are just comparing all the quads looking for the
    # one that has the best score. This has to be at least 10.0 as, based on
    # testing, you can sometimes get scores as high as X.X even with two random
    # data sets.
    #
    best_ratio = 10.0
    best_transform = None
    matches = 0
    for q1 in quads1:
        for q2 in quads2:
            if q1.isMatch(q2, tolerance = args.tolerance):
                fg_p = fgProbability(kd1, kd2, q1.getTransform(q2))
                ratio = math.log(fg_p/bg_p)
                print("Match", matches, fg_p, bg_p, ratio)
                if (ratio > best_ratio):
                    best_ratio = ratio
                    best_transform = q1.getTransform(q2) + q2.getTransform(q1)
                matches += 1

    print("Found", matches, "matching quads")

    if best_transform is not None:
        plotMatch(kd1, kd2, best_transform, args.results + ".png", show = (not args.no_plots))

        #
        # Save mapping using the same format that multi-plane uses.
        #
        mapping = {"1_0_x" : best_transform[0],
                   "1_0_y" : best_transform[1],
                   "0_1_x" : best_transform[2],
                   "0_1_y" : best_transform[3]}

        with open(args.results, 'wb') as fp:
            pickle.dump(mapping, fp)

    else:
        print("No transform of sufficient quality was found.")

