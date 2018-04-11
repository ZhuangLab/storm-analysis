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
images, The Astronomical Journal 139, 1782-1800.

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
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.micrometry.quads as quads


def applyTransform(kd, transform):
    tx = transform[0]
    ty = transform[1]
    x = tx[0] + tx[1]*kd.data[:,0] + tx[2]*kd.data[:,1]
    y = ty[0] + ty[1]*kd.data[:,0] + ty[2]*kd.data[:,1]
    return [x, y]


def fgProbability(kd1, kd2, transform, bg_p):
    """
    Returns an estimate of how likely the transform is correct.
    """
    # Transform 'other' coordinates into the 'reference' frame.
    [x2, y2] = applyTransform(kd2, transform)
    p2 = numpy.stack((x2, y2), axis = -1)

    # Calculate distance to nearest point in 'reference'.
    [dist, index] = kd1.query(p2)

    # Score assuming a localization accuracy of 1 pixel.
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


def makeTreeAndQuadsFromH5File(h5_filename, min_size = None, max_size = None, max_neighbors = 10):
    """
    Make a KD tree and a list of quads from a SA HDF5 file.

    Note: This file should probably only have localizations for a single frame.
    """
    with saH5Py.SAH5Py(h5_filename) as h5:
        locs = h5.getLocalizations()

        # Check that there is only a single frame with localizations.
        cnts = 0
        for elt in h5.localizationsIterator():
            cnts += 1

        # Warning if there is more than 1 frame in the data.
        if (cnts > 1):
            print("Warning: Localizations in multiple frames detected!")

    return makeTreeAndQuads(locs['x'],
                            locs['y'],
                            min_size = min_size,
                            max_size = max_size,
                            max_neighbors = max_neighbors)


def plotMatch(kd1, kd2, transform, save_as = None, show = True):
    [x2, y2] = applyTransform(kd2, transform)
    
    fig = pyplot.figure()
    pyplot.scatter(kd1.data[:,0], kd1.data[:,1], facecolors = 'none', edgecolors = 'red', s = 100)
    pyplot.scatter(x2, y2, color = 'green', marker = '+', s = 100)

    legend = pyplot.legend(('reference', 'other'), loc=1)
    pyplot.xlabel("pixels")
    pyplot.ylabel("pixels")

    ax = pyplot.gca()
    ax.set_aspect('equal')

    if save_as is not None:
        fig.savefig(save_as)
    
    if show:
        pyplot.show()

def prettyPrintTransform(transform):
    """
    Pretty print the transform.
    """
    print("Reference to other transform:")
    print("  {0:.4f} {1:.4f} {2:.4f}".format(transform[2][0], transform[2][1], transform[2][2]))
    print("  {0:.4f} {1:.4f} {2:.4f}".format(transform[3][0], transform[3][1], transform[3][2]))
    print("")
    print("Other to reference transform:")
    print("  {0:.4f} {1:.4f} {2:.4f}".format(transform[0][0], transform[0][1], transform[0][2]))
    print("  {0:.4f} {1:.4f} {2:.4f}".format(transform[1][0], transform[1][1], transform[1][2]))
    print("")
    
def runMicrometry(locs1, locs2, results, min_size = 5.0, max_size = 100.0, max_neighbors = 20, tolerance = 1.0e-2, no_plots = False):
    """
    This performs all the steps in micrometry in a single function.

    locs1 - The name of the first HDF5 format file containing the localizations.
    locs2 - The name of the second HDF5 format file containing the localizations.
    results - The name of the file to save the mapping in.
    min_size - The minimum quad size in pixels.
    max_size - The maximum quad size in pixels.
    max_neighbors - The maximum number of (nearest) neighbors to consider in quad formation.
    tolerance - Tolerance is a relative unit used in quad matching.
    """
    mm = Micrometry(locs1,
                    min_size = min_size,
                    max_size = max_size,
                    max_neighbors = max_neighbors)
    [best_ratio, best_transform] = mm.findTransform(locs2, tolerance)

    if (best_ratio > 10.0):
        plotMatch(mm.getRefKDTree(),
                  mm.getOtherKDTree(),
                  best_transform,
                  save_as = results + ".png",
                  show = (not no_plots))

        #
        # Save mapping using the same format that multi-plane uses.
        #
        mapping = {"1_0_x" : best_transform[0],
                   "1_0_y" : best_transform[1],
                   "0_1_x" : best_transform[2],
                   "0_1_y" : best_transform[3]}

        with open(results, 'wb') as fp:
            pickle.dump(mapping, fp)

    else:
        print("No transform of sufficient quality was found.")
        if best_transform is not None:
            plotMatch(mm.getRefKDTree(),
                      mm.getOtherKDTree(),
                      best_transform,
                      show = (not no_plots))


class Micrometry(object):
    """
    Class for performing geometric hashing to identify the affine 
    transform between two localization files.
    """
    def __init__(self, ref_filename = None, min_size = None, max_size = None, max_neighbors = None, verbose = True, **kwds):
        """
        ref_filename - The name of an HDF5 format file containing the localizations.
        min_size - The minimum quad size in pixels.
        max_size - The maximum quad size in pixels.
        max_neighbors - The maximum number of (nearest) neighbors to consider in quad formation.
        """
        super(Micrometry, self).__init__(**kwds)

        self.density = None
        self.kd_other = None
        self.kd_ref = None
        self.max_neighbors = max_neighbors
        self.max_size = max_size
        self.min_size = min_size
        self.quads_other = None
        self.quads_ref = None
        self.verbose = verbose

        # Create quads for the reference data.
        #
        if self.verbose:
            print("Making quads for the 'reference' data.")
        [self.kd_ref, self.quads_ref] = makeTreeAndQuadsFromH5File(ref_filename,
                                                                   min_size = self.min_size,
                                                                   max_size = self.max_size,
                                                                   max_neighbors = self.max_neighbors)
        if self.verbose:
            print("Created", len(self.quads_ref), "quads")
            print("")

        # Estimate background, the density of points in the reference.
        #
        with saH5Py.SAH5Py(ref_filename) as h5:
            [movie_x, movie_y] = h5.getMovieInformation()[:2]
        self.density = 1.0/(movie_x * movie_y)

    def getOtherKDTree(self):
        return self.kd_other
    
    def getRefKDTree(self):
        return self.kd_ref

    def findTransform(self, other_filename, tolerance, min_size = None, max_size = None, max_neighbors = None):

        if max_neighbors is None:
            max_neighbors = self.max_neighbors
            
        if max_size is None:
            max_size = self.max_size
            
        if min_size is None:
            min_size = self.min_size

        # Create quads for the other data.
        #
        if self.verbose:
            print("Making quads for the 'other' data.")
        [self.kd_other, self.quads_other] = makeTreeAndQuadsFromH5File(other_filename,
                                                                       min_size = min_size,
                                                                       max_size = max_size,
                                                                       max_neighbors = max_neighbors)

        if self.verbose:
            print("Created", len(self.quads_other), "quads")
            print("")
            print("Comparing quads.")
        
        #
        # Unlike astrometry.net we are just comparing all the quads looking for the
        # one that has the best score. This should be at least 10.0 as, based on
        # testing, you can sometimes get scores as high as 9.7 even if the match
        # is not actually any good.
        #
        best_ratio = 0.0
        best_transform = None
        matches = 0
        for q1 in self.quads_ref:
            for q2 in self.quads_other:
                if q1.isMatch(q2, tolerance = tolerance):
                    fg_p = fgProbability(self.kd_ref, self.kd_other, q1.getTransform(q2), self.density)
                    ratio = math.log(fg_p/self.density)
                    if self.verbose:
                        print("Match {0:d} {1:.2f} {2:.2e} {3:.2f}".format(matches, fg_p, self.density, ratio))
                    if (ratio > best_ratio):
                        best_ratio = ratio
                        best_transform = q1.getTransform(q2) + q2.getTransform(q1)
                    matches += 1

        if self.verbose:
            print("Found", matches, "matching quads")

        return [best_ratio, best_transform]
    
    
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
                        help = "Maximum neighbors to search when making quads, default is 20")
    parser.add_argument('--tolerance', dest='tolerance', type=float, required=False, default=1.0e-2,
                        help = "Tolerance for matching quads, default is 1.0e-2.")
    parser.add_argument('--no_plots', dest='no_plots', action='store_true', default=False,
                        help = "Don't show a plot of the results.")

    args = parser.parse_args()

    runMicrometry(args.locs1,
                  args.locs2,
                  args.results,                  
                  min_size = args.min_size,
                  max_size = args.max_size,
                  max_neighbors = args.max_neighbors,
                  tolerance = args.tolerance,
                  no_plots = args.no_plots)
