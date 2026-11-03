#!/usr/bin/env python
"""
Given two lists of localizations, returns the transform between
them. This is degree 1 affine transform. As such this could fail
for images with large differences in field curvature.

This uses the ideas in this paper, applied to fiducial references
like flourescent beads:

Lang, D., Hogg, D. W., Mierle, K., Blanton, M., & Roweis, S., 2010,
Astrometry.net: Blind astrometric calibration of arbitrary astronomical
images, The Astronomical Journal 139, 1782-1800.

The steps are:

1. Make quads for both lists. The 'other' list is also used with
   x negated, so that mirrored images match.
2. Pair 'reference' and 'other' quads with similar codes, most
   similar first.
3. Fit an affine transform to the four points of each pair, and
   count how many 'other' localizations it puts within the pairing
   radius of a 'reference' localization.
4. Refine the best transforms by least squares on all the pairs.
5. Score the best transform by how unlikely its number of pairs
   would be if the two lists were unrelated.

The pairing radius is max(1, 3 * sqrt(2) * localization_error)
pixels, enough for the pair separation of two localizations that
each have this error in x and y.

The score is -ln of the chance of finding this many pairs or more
with a random transform, multiplied by the number of transforms
that were tried. Scores above 10 are significant whatever the size
of the field or the number of localizations.

Hazen 07/17
"""
import math
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import pickle
import scipy
import scipy.spatial
import scipy.stats

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.micrometry.quads as quads


def applyTransform(kd, transform):
    tx = transform[0]
    ty = transform[1]
    x = tx[0] + tx[1]*kd.data[:,0] + tx[2]*kd.data[:,1]
    y = ty[0] + ty[1]*kd.data[:,0] + ty[2]*kd.data[:,1]
    return [x, y]


def boundingBox(points):
    """
    Returns [x_min, y_min, x_max, y_max] for a (N, 2) array of points.
    """
    return [numpy.min(points[:,0]), numpy.min(points[:,1]), numpy.max(points[:,0]), numpy.max(points[:,1])]


def chancePairs(kd_ref, other_points, transform, radius):
    """
    Returns the number of pairs within radius to expect if the transform
    were wrong. This is the number of transformed 'other' localizations
    close enough to the 'reference' field to pair, times the number of
    'reference' localizations to expect within radius of each, given
    their density where the two fields overlap.

    The fields are the bounding boxes of the localizations. The bounding
    box is used rather than the size of the movie as it is never larger
    than the area the localizations actually cover, so the estimate of
    the number of chance pairs errs on the high side.
    """
    ref_box = boundingBox(kd_ref.data)
    other_box = boundingBox(other_points)

    # Overlap of the 'reference' field and the transformed 'other' field.
    corners = numpy.array([[other_box[0], other_box[1]],
                           [other_box[2], other_box[1]],
                           [other_box[2], other_box[3]],
                           [other_box[0], other_box[3]]])
    overlap = clipToBox(transformPoints(corners, transform), ref_box)
    area = polygonArea(overlap)
    if not (area > 1.0e-6):
        return numpy.inf

    # Density of 'reference' localizations in the overlap.
    ref_in = transformPoints(kd_ref.data, invertTransform(transform))
    ref_density = numpy.count_nonzero(insideBox(ref_in, other_box)) / area

    # 'Other' localizations close enough to the 'reference' field to pair.
    near_box = [ref_box[0] - radius, ref_box[1] - radius, ref_box[2] + radius, ref_box[3] + radius]
    n_other = numpy.count_nonzero(insideBox(transformPoints(other_points, transform), near_box))

    return n_other * ref_density * math.pi * radius * radius


def clipToBox(polygon, box):
    """
    Clip a convex polygon to a box [x_min, y_min, x_max, y_max]
    (Sutherland-Hodgman). Returns the vertices of the clipped polygon.
    """
    polygon = list(polygon)
    for [axis, limit, sign] in [[0, box[0], 1.0], [0, box[2], -1.0], [1, box[1], 1.0], [1, box[3], -1.0]]:
        clipped = []
        for i in range(len(polygon)):
            current = polygon[i]
            previous = polygon[i-1]
            current_inside = (sign * (current[axis] - limit) >= 0.0)
            previous_inside = (sign * (previous[axis] - limit) >= 0.0)
            if (current_inside != previous_inside):
                t = (limit - previous[axis])/(current[axis] - previous[axis])
                clipped.append(previous + t * (current - previous))
            if current_inside:
                clipped.append(current)
        polygon = clipped
        if (len(polygon) == 0):
            break
    return polygon


def countDistinctPartners(kd_ref, probe_points, transforms, radius):
    """
    For each of the (C, 2, 3) transforms, count the 'reference'
    localizations that are the nearest neighbor of a transformed
    probe point within radius. Each 'reference' localization is
    counted once however many probe points it is the nearest to.
    """
    px = probe_points[:,0][None,:]
    py = probe_points[:,1][None,:]
    x = transforms[:,0,0][:,None] + transforms[:,0,1][:,None]*px + transforms[:,0,2][:,None]*py
    y = transforms[:,1,0][:,None] + transforms[:,1,1][:,None]*px + transforms[:,1,2][:,None]*py

    n_ref = kd_ref.data.shape[0]
    [dist, index] = kd_ref.query(numpy.stack((x.ravel(), y.ravel()), axis = -1), distance_upper_bound = radius)

    # Misses have index n_ref.
    index = numpy.sort(index.reshape(x.shape), axis = 1)
    is_new = numpy.ones(index.shape, dtype = bool)
    is_new[:,1:] = (index[:,1:] != index[:,:-1])
    return numpy.count_nonzero(is_new & (index < n_ref), axis = 1)


def findPairs(kd_ref, other_points, transform, radius):
    """
    Pair transformed 'other' localizations with 'reference' localizations
    within radius, one to one, closest pairs first.

    Returns [ref_index, other_index].
    """
    n_ref = kd_ref.data.shape[0]
    k = min(3, n_ref)
    [dist, index] = kd_ref.query(transformPoints(other_points, transform), k = k, distance_upper_bound = radius)
    dist = dist.reshape(other_points.shape[0], k)
    index = index.reshape(other_points.shape[0], k)

    mask = numpy.isfinite(dist)
    other_index = numpy.repeat(numpy.arange(other_points.shape[0])[:,None], k, axis = 1)[mask]
    ref_index = index[mask]
    order = numpy.argsort(dist[mask], kind = "stable")

    ref_used = numpy.zeros(n_ref, dtype = bool)
    other_used = numpy.zeros(other_points.shape[0], dtype = bool)
    ref_paired = []
    other_paired = []
    for i in order:
        if ref_used[ref_index[i]] or other_used[other_index[i]]:
            continue
        ref_used[ref_index[i]] = True
        other_used[other_index[i]] = True
        ref_paired.append(ref_index[i])
        other_paired.append(other_index[i])

    return [numpy.array(ref_paired, dtype = numpy.int64), numpy.array(other_paired, dtype = numpy.int64)]


def findTransformFromPoints(ref_points, other_points, min_size = None, max_size = None, max_neighbors = 8, tolerance = 4.0e-2, localization_error = 0.5, candidate_budget = 2000, verbose = False):
    """
    Find the transform between two (N, 2) arrays of x, y positions.

    See Micrometry.findTransform() for the parameters and the results.
    """
    mm = Micrometry(ref_points = ref_points,
                    min_size = min_size,
                    max_size = max_size,
                    max_neighbors = max_neighbors,
                    verbose = verbose)
    return mm.matchPoints(other_points,
                          tolerance = tolerance,
                          localization_error = localization_error,
                          candidate_budget = candidate_budget)


def fitAffine(source, target):
    """
    Returns the affine transform, as [tx, ty], that maps the (N, 2)
    source points onto the (N, 2) target points by least squares.
    """
    m = numpy.ones((source.shape[0], 3))
    m[:,1:] = source
    params = numpy.linalg.lstsq(m, target, rcond = None)[0]
    return [params[:,0], params[:,1]]


def fitAffines(sources, targets):
    """
    Returns the (C, 2, 3) affine transforms that map each set of (C, 4, 2)
    source points onto the matching target points by least squares.
    """
    m = numpy.ones(sources.shape[:2] + (3,))
    m[:,:,1:] = sources
    params = numpy.matmul(numpy.linalg.pinv(m), targets)
    return numpy.transpose(params, (0, 2, 1))


def insideBox(points, box):
    return ((points[:,0] >= box[0]) & (points[:,0] <= box[2]) & (points[:,1] >= box[1]) & (points[:,1] <= box[3]))


def invertTransform(transform):
    """
    Returns the inverse of the affine transform [tx, ty].
    """
    matrix = numpy.array([[transform[0][1], transform[0][2], transform[0][0]],
                          [transform[1][1], transform[1][2], transform[1][0]],
                          [0.0, 0.0, 1.0]])
    inverse = numpy.linalg.inv(matrix)
    return [inverse[0,[2,0,1]], inverse[1,[2,0,1]]]


def loadLocalizations(h5_filename):
    """
    Returns the localizations in a SA HDF5 file as a (N, 2) array.

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

    return numpy.stack((locs['x'], locs['y']), axis = -1).astype(numpy.float64)


def pairingRadius(localization_error):
    """
    The pairing radius in pixels for localizations that each have
    localization_error pixels of error in x and y.
    """
    return max(1.0, 3.0 * math.sqrt(2.0) * localization_error)


def plotMatch(kd1, kd2, transform, save_as = None, show = True):
    [x2, y2] = applyTransform(kd2, transform)

    fig = pyplot.figure(figsize = (10,10))
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


def polygonArea(polygon):
    if (len(polygon) < 3):
        return 0.0
    polygon = numpy.array(polygon)
    x = polygon[:,0]
    y = polygon[:,1]
    return 0.5 * abs(numpy.dot(x, numpy.roll(y, -1)) - numpy.dot(y, numpy.roll(x, -1)))


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


def refineTransform(kd_ref, other_points, transform, radius, max_rounds = 10):
    """
    Refine the 'other' to 'reference' transform by least squares on the
    one to one pairs within radius. This repeats until the pairs stop
    changing, or for max_rounds.

    Returns [transform, [ref_index, other_index]].
    """
    pairs = findPairs(kd_ref, other_points, transform, radius)
    for i in range(max_rounds):
        if (pairs[0].size < 4):
            break
        transform = fitAffine(other_points[pairs[1]], kd_ref.data[pairs[0]])
        new_pairs = findPairs(kd_ref, other_points, transform, radius)
        if (new_pairs[0].size == pairs[0].size) and numpy.array_equal(numpy.sort(new_pairs[1]), numpy.sort(pairs[1])):
            pairs = new_pairs
            break
        pairs = new_pairs

    return [transform, pairs]


def reverseTransform(kd_ref, other_points, transform, pairs):
    """
    Returns the 'reference' to 'other' transform, fit by least squares
    to the pairs from refineTransform() if there are enough of them,
    otherwise the inverse of transform.
    """
    if (pairs[0].size >= 4):
        return fitAffine(kd_ref.data[pairs[0]], other_points[pairs[1]])
    else:
        return invertTransform(transform)


def runMicrometry(locs1, locs2, results, min_size = None, max_size = None, max_neighbors = 8, tolerance = 4.0e-2, no_plots = False, localization_error = 0.5, candidate_budget = 2000):
    """
    This performs all the steps in micrometry in a single function.

    locs1 - The name of the first HDF5 format file containing the localizations.
    locs2 - The name of the second HDF5 format file containing the localizations.
    results - The name of the file to save the mapping in.
    min_size - The minimum quad size in pixels, default is no limit.
    max_size - The maximum quad size in pixels, default is no limit.
    max_neighbors - The number of (nearest) neighbors to group with each localization in quad formation.
    tolerance - Tolerance is a relative unit used in quad matching.
    no_plots - Don't show a plot of the results.
    localization_error - The localization error in pixels.
    candidate_budget - The maximum number of quad matches to check.

    Returns [score, transform].
    """
    mm = Micrometry(locs1,
                    min_size = min_size,
                    max_size = max_size,
                    max_neighbors = max_neighbors)
    [score, transform] = mm.findTransform(locs2,
                                          tolerance,
                                          localization_error = localization_error,
                                          candidate_budget = candidate_budget)

    if (score > 10.0):
        plotMatch(mm.getRefKDTree(),
                  mm.getOtherKDTree(),
                  transform,
                  save_as = results + ".png",
                  show = (not no_plots))

        saveMapping(transform, results)

    else:
        print("No transform of sufficient quality was found.")
        if transform is not None:
            plotMatch(mm.getRefKDTree(),
                      mm.getOtherKDTree(),
                      transform,
                      show = (not no_plots))

    return [score, transform]


def saveMapping(transform, filename):
    """
    Save a transform from findTransform() in the mapping format
    that multi-plane uses. '0' is the reference, '1' is the other.
    """
    mapping = {"1_0_x" : transform[0],
               "1_0_y" : transform[1],
               "0_1_x" : transform[2],
               "0_1_y" : transform[3]}

    with open(filename, 'wb') as fp:
        pickle.dump(mapping, fp)


def significance(n_pairs, chance, n_tests):
    """
    Returns -ln of the probability of n_pairs or more pairs when chance
    pairs are expected, multiplied by the number of tests, or 0.0 if
    that probability is more than 1.
    """
    if (n_pairs <= 0) or not numpy.isfinite(chance):
        return 0.0
    chance = max(chance, 1.0e-12)

    log_tail = scipy.stats.poisson.logsf(n_pairs - 1, chance)
    if not numpy.isfinite(log_tail):
        # The tail is too small for logsf. For n_pairs well above
        # chance it is pmf(n_pairs) * (1 + chance/(n_pairs + 1) + ...),
        # which is at most this.
        log_tail = scipy.stats.poisson.logpmf(n_pairs, chance) - math.log1p(-chance/(n_pairs + 1))

    return max(0.0, -(log_tail + math.log(n_tests)))


def transformPoints(points, transform):
    """
    Apply the affine transform [tx, ty] to a (N, 2) array of points.
    """
    tx = transform[0]
    ty = transform[1]
    return numpy.stack((tx[0] + tx[1]*points[:,0] + tx[2]*points[:,1],
                        ty[0] + ty[1]*points[:,0] + ty[2]*points[:,1]), axis = -1)


class Micrometry(object):
    """
    Class for performing geometric hashing to identify the affine
    transform between two localization files.
    """
    def __init__(self, ref_filename = None, min_size = None, max_size = None, max_neighbors = 8, verbose = True, ref_points = None, **kwds):
        """
        ref_filename - The name of an HDF5 format file containing the localizations.
        min_size - The minimum quad size in pixels, default is no limit.
        max_size - The maximum quad size in pixels, default is no limit.
        max_neighbors - The number of (nearest) neighbors to group with each localization in quad formation.
        verbose - Print progress.
        ref_points - A (N, 2) array of x, y positions to use instead of ref_filename.
        """
        super(Micrometry, self).__init__(**kwds)

        self.kd_other = None
        self.kd_ref = None
        self.max_neighbors = max_neighbors
        self.max_size = max_size
        self.min_size = min_size
        self.quads_ref = None
        self.verbose = verbose

        if ref_points is None:
            ref_points = loadLocalizations(ref_filename)

        # Create quads for the reference data.
        #
        if self.verbose:
            print("Making quads for the 'reference' data.")
        self.kd_ref = scipy.spatial.KDTree(numpy.asarray(ref_points, dtype = numpy.float64))
        self.quads_ref = quads.makeQuads(self.kd_ref.data,
                                         min_size = self.min_size,
                                         max_size = self.max_size,
                                         max_neighbors = self.max_neighbors)
        if self.verbose:
            print("Created", self.quads_ref[0].shape[0], "quads")
            print("")

    def getOtherKDTree(self):
        return self.kd_other

    def getRefKDTree(self):
        return self.kd_ref

    def findTransform(self, other_filename, tolerance = 4.0e-2, min_size = None, max_size = None, max_neighbors = None, localization_error = 0.5, candidate_budget = 2000):
        """
        Find the transform between the 'reference' localizations and
        those in other_filename.

        other_filename - The name of an HDF5 format file containing the localizations.
        tolerance - Quads match if no element of their codes differs by more than this.
        min_size, max_size, max_neighbors - Quad formation for the 'other' localizations,
                     default is the same as for the 'reference' localizations.
        localization_error - The localization error in pixels.
        candidate_budget - The maximum number of quad matches to check, closest codes first.

        Returns [score, transform]. A score above 10 is significant (see the
        module docstring). transform is [other to reference x, other to
        reference y, reference to other x, reference to other y], each as
        [c0, cx, cy] so that x' = c0 + cx * x + cy * y, or None if no quads
        matched.
        """
        return self.matchPoints(loadLocalizations(other_filename),
                                tolerance = tolerance,
                                min_size = min_size,
                                max_size = max_size,
                                max_neighbors = max_neighbors,
                                localization_error = localization_error,
                                candidate_budget = candidate_budget)

    def matchPoints(self, other_points, tolerance = 4.0e-2, min_size = None, max_size = None, max_neighbors = None, localization_error = 0.5, candidate_budget = 2000):
        """
        Find the transform between the 'reference' localizations and
        a (N, 2) array of x, y positions. See findTransform().
        """
        if max_neighbors is None:
            max_neighbors = self.max_neighbors

        if max_size is None:
            max_size = self.max_size

        if min_size is None:
            min_size = self.min_size

        other_points = numpy.asarray(other_points, dtype = numpy.float64)
        self.kd_other = scipy.spatial.KDTree(other_points)
        radius = pairingRadius(localization_error)

        # Create quads for the other data, as is and mirrored.
        #
        if self.verbose:
            print("Making quads for the 'other' data.")
        [other_indices, other_codes] = quads.makeQuads(other_points,
                                                       min_size = min_size,
                                                       max_size = max_size,
                                                       max_neighbors = max_neighbors)
        [mirror_indices, mirror_codes] = quads.makeQuads(other_points * numpy.array([-1.0, 1.0]),
                                                         min_size = min_size,
                                                         max_size = max_size,
                                                         max_neighbors = max_neighbors)
        other_indices = numpy.concatenate((other_indices, mirror_indices))
        other_codes = numpy.concatenate((other_codes, mirror_codes))

        if self.verbose:
            print("Created", other_indices.shape[0], "quads, including mirrored quads")
            print("")
            print("Comparing quads.")

        [ref_indices, ref_codes] = self.quads_ref
        if (ref_indices.shape[0] == 0) or (other_indices.shape[0] == 0):
            if self.verbose:
                print("Not enough quads to compare.")
            return [0.0, None]

        #
        # Quad pairs whose codes differ by at most tolerance in each
        # element, closest codes first.
        #
        matches = scipy.spatial.cKDTree(ref_codes).sparse_distance_matrix(scipy.spatial.cKDTree(other_codes),
                                                                          tolerance,
                                                                          p = numpy.inf,
                                                                          output_type = "ndarray")
        if (matches.size == 0):
            if self.verbose:
                print("No matching quads.")
            return [0.0, None]

        code_distance = numpy.linalg.norm(ref_codes[matches["i"]] - other_codes[matches["j"]], axis = 1)
        order = numpy.argsort(code_distance, kind = "stable")[:candidate_budget]
        ref_quads = ref_indices[matches["i"][order]]
        other_quads = other_indices[matches["j"][order]]

        #
        # Fit an affine transform to each pair of quads. The mirrored
        # quads use the original positions, so the fit includes the
        # mirroring.
        #
        transforms = fitAffines(other_points[other_quads], self.kd_ref.data[ref_quads])
        good = numpy.all(numpy.isfinite(transforms), axis = (1,2))
        transforms = transforms[good]
        n_tests = transforms.shape[0]
        if (n_tests == 0):
            return [0.0, None]

        #
        # Rank the transforms by how many localizations they pair,
        # using a subset of the 'other' localizations for speed.
        #
        rng = numpy.random.default_rng(0)
        probe_size = min(50, other_points.shape[0])
        probe_points = other_points[rng.choice(other_points.shape[0], size = probe_size, replace = False)]
        counts = countDistinctPartners(self.kd_ref, probe_points, transforms, radius)
        top = numpy.argsort(-counts, kind = "stable")[:5]

        #
        # Refine the top transforms using all the localizations, and
        # keep the most significant.
        #
        best = None
        for i in top:
            transform = [transforms[i,0,:], transforms[i,1,:]]
            [transform, pairs] = refineTransform(self.kd_ref, other_points, transform, radius)

            # The four points of the quad pair are paired whether or not the
            # transform is right, so they are not counted as evidence.
            n_pairs = pairs[0].size
            score = significance(n_pairs - 4, chancePairs(self.kd_ref, other_points, transform, radius), n_tests)
            if best is None or ([score, n_pairs] > [best[0], best[1]]):
                best = [score, n_pairs, transform, pairs]

        [score, n_pairs, transform, pairs] = best
        reverse = reverseTransform(self.kd_ref, other_points, transform, pairs)

        if self.verbose:
            print("Checked", n_tests, "quad matches")
            print("Best match pairs {0:d} localizations, score {1:.1f}".format(n_pairs, score))
            if ((transform[0][1] * transform[1][2] - transform[0][2] * transform[1][1]) < 0.0):
                print("The 'other' localizations are mirrored")
            print("")

        return [score, list(transform) + list(reverse)]


if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = 'Micrometry - ...')

    parser.add_argument('--locs1', dest='locs1', type=str, required=True,
                        help = "The name of the 'reference' localizations file")
    parser.add_argument('--locs2', dest='locs2', type=str, required=True,
                        help = "The name of the 'other' localizations file")
    parser.add_argument('--results', dest='results', type=str, required=True,
                        help = "The name of the file to save the transform (if any) in.")
    parser.add_argument('--min_size', dest='min_size', type=float, required=False, default=None,
                        help = "Minimum quad size (pixels), default is no limit.")
    parser.add_argument('--max_size', dest='max_size', type=float, required=False, default=None,
                        help = "Maximum quad size (pixels), default is no limit.")
    parser.add_argument('--max_neighbors', dest='max_neighbors', type=int, required=False, default=8,
                        help = "Number of nearest neighbors to group with each localization when making quads, default is 8")
    parser.add_argument('--tolerance', dest='tolerance', type=float, required=False, default=4.0e-2,
                        help = "Tolerance for matching quads, default is 4.0e-2.")
    parser.add_argument('--localization_error', dest='localization_error', type=float, required=False, default=0.5,
                        help = "Localization error (pixels) in x and y, default is 0.5.")
    parser.add_argument('--candidate_budget', dest='candidate_budget', type=int, required=False, default=2000,
                        help = "Maximum number of quad matches to check, default is 2000.")
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
                  no_plots = args.no_plots,
                  localization_error = args.localization_error,
                  candidate_budget = args.candidate_budget)
