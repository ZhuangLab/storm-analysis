#!/usr/bin/env python
"""
Quads, groups of four points described by a code that does not
change under translation, rotation or scaling.

A quad's code is the position of its two inner points C and D in a
frame where A is at (0, 0) and B is at (1, 1). A and B are the two
points furthest apart, and C and D must lie inside the circle with
diameter AB. Swapping A and B, or C and D, describes the same four
points differently, so a quad is labelled such that xc <= xd and
xc + xd <= 1. That gives each set of four points a single code.

A mirror image of a quad does not have the same code. To match
mirrored point sets, make the quads of one of them a second time
with x negated (see micrometry.py).

Hazen 07/17
"""
import itertools
import math
import numpy
import scipy
import scipy.spatial

c45 = math.cos(0.25 * math.pi)
s45 = math.sin(0.25 * math.pi)
root2 = math.sqrt(2.0)

# The pairs of points in a group of four, and the other two points for each pair.
pairs = numpy.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]])
pair_complements = numpy.array([[2,3], [1,3], [1,2], [0,3], [0,2], [0,1]])


def makeQuad(A, B, C, D):
    """
    Returns the code [xc, yc, xd, yd] if points A,B,C,D form a
    proper quad with this labelling, otherwise returns None.

    A,B define the coordinate system of the quad.
    C,D are the internal points.

    This is the reference definition of the code, makeQuads()
    computes the same thing for many quads at once.
    """

    # Calculate scale.
    dab_x = B[0] - A[0]
    dab_y = B[1] - A[1]
    dab_l = math.sqrt(dab_x*dab_x + dab_y*dab_y)

    dab_l = 1.0/dab_l

    # Calculate circle center.
    cx = 0.5*(A[0]+B[0])
    cy = 0.5*(A[1]+B[1])

    # Calculate radius (squared).
    dx = A[0] - cx
    dy = A[1] - cy
    max_rr = dx*dx + dy*dy

    # Verify that C,B are within radius of the center point.
    for P in [C,D]:
        dx = P[0] - cx
        dy = P[1] - cy
        rr = dx*dx + dy*dy
        if (rr > max_rr):
            return

    # Calculate basis vectors.
    dab_x = dab_x * dab_l
    dab_y = dab_y * dab_l

    x_vec = [c45 * dab_x + s45 * dab_y, -s45 * dab_x + c45*dab_y]
    y_vec = [c45 * dab_x - s45 * dab_y, s45 * dab_x + c45*dab_y]

    dab_l = root2 * dab_l

    # Calculate xc, yc.
    dac_x = dab_l * (C[0] - A[0])
    dac_y = dab_l * (C[1] - A[1])

    xc = x_vec[0] * dac_x + x_vec[1] * dac_y
    yc = y_vec[0] * dac_x + y_vec[1] * dac_y

    # Calcule xd, yd.
    dad_x = dab_l * (D[0] - A[0])
    dad_y = dab_l * (D[1] - A[1])

    xd = x_vec[0] * dad_x + x_vec[1] * dad_y
    yd = y_vec[0] * dad_x + y_vec[1] * dad_y

    if (xc > xd):
        return

    if ((xc + xd) > 1.0):
        return

    return [xc, yc, xd, yd]


def makeQuads(points, min_size = None, max_size = None, max_neighbors = 8):
    """
    Construct the quads for a set of points.

    Every group of four drawn from a point and its max_neighbors
    nearest neighbors is a candidate quad. A group appears once
    however many neighborhoods it is part of.

    Note: The number of candidate quads grows as the number of
          points times max_neighbors to the 4th power.

    points - A (N, 2) array of x, y positions.
    min_size - A,B points must be at least this distance from each
               other, default is no limit.
    max_size - A,B points must be at most this distance from each
               other, default is no limit.
    max_neighbors - The number of nearest neighbors to group with
               each point, default is 8.

    Returns [indices, codes]. indices is a (M, 4) array with the
    index of each quad's A, B, C and D points in points, codes is
    a (M, 4) array with each quad's xc, yc, xd, yd.
    """
    points = numpy.asarray(points, dtype = numpy.float64)
    if (points.shape[0] < 4):
        return [numpy.zeros((0, 4), dtype = numpy.int64), numpy.zeros((0, 4))]

    #
    # Each point and its nearest neighbors. A point is always its
    # own nearest neighbor.
    #
    group_size = min(max_neighbors + 1, points.shape[0])
    kd = scipy.spatial.cKDTree(points)
    [dist, neighborhoods] = kd.query(points, k = group_size)

    #
    # All the groups of four in each neighborhood, without repeats.
    #
    choices = numpy.array(list(itertools.combinations(range(group_size), 4)))
    groups = neighborhoods[:, choices].reshape(-1, 4)
    groups = numpy.unique(numpy.sort(groups, axis = 1), axis = 0)
    group_points = points[groups]

    #
    # A and B are the pair furthest apart, C and D are the other two.
    #
    deltas = group_points[:, pairs[:,0], :] - group_points[:, pairs[:,1], :]
    pair_rr = numpy.sum(deltas * deltas, axis = 2)
    widest = numpy.argmax(pair_rr, axis = 1)
    rows = numpy.arange(groups.shape[0])
    labels = numpy.concatenate((pairs[widest], pair_complements[widest]), axis = 1)
    indices = groups[rows[:, None], labels]
    ab_rr = pair_rr[rows, widest]

    #
    # Keep groups whose A,B separation is within the size limits.
    #
    mask = (ab_rr > 0.0)
    if min_size is not None:
        mask = mask & (ab_rr >= min_size * min_size)
    if max_size is not None:
        mask = mask & (ab_rr <= max_size * max_size)

    #
    # Keep groups where C and D are inside the circle with diameter AB.
    #
    center = 0.5 * (points[indices[:,0]] + points[indices[:,1]])
    for k in [2, 3]:
        dcenter = points[indices[:,k]] - center
        mask = mask & (numpy.sum(dcenter * dcenter, axis = 1) <= 0.25 * ab_rr)

    indices = indices[mask]
    ab_rr = ab_rr[mask]

    #
    # Code of C and D in the frame where A is (0, 0) and B is (1, 1).
    #
    A = points[indices[:,0]]
    dab = points[indices[:,1]] - A
    x_axis = numpy.stack((dab[:,0] + dab[:,1], dab[:,1] - dab[:,0]), axis = 1) / ab_rr[:, None]
    y_axis = numpy.stack((dab[:,0] - dab[:,1], dab[:,0] + dab[:,1]), axis = 1) / ab_rr[:, None]

    dac = points[indices[:,2]] - A
    dad = points[indices[:,3]] - A
    codes = numpy.stack((numpy.sum(x_axis * dac, axis = 1),
                         numpy.sum(y_axis * dac, axis = 1),
                         numpy.sum(x_axis * dad, axis = 1),
                         numpy.sum(y_axis * dad, axis = 1)), axis = 1)

    #
    # Label so that xc + xd <= 1. Swapping A and B takes each code
    # point (x, y) to (1 - x, 1 - y).
    #
    swap = ((codes[:,0] + codes[:,2]) > 1.0)
    codes[swap] = 1.0 - codes[swap]
    indices[swap] = indices[swap][:, [1, 0, 2, 3]]

    #
    # Label so that xc <= xd.
    #
    swap = (codes[:,0] > codes[:,2])
    codes[swap] = codes[swap][:, [2, 3, 0, 1]]
    indices[swap] = indices[swap][:, [0, 1, 3, 2]]

    return [indices, codes]


if (__name__ == "__main__"):

    import time

    numpy.random.seed(0)
    xp = numpy.random.uniform(low = 0.0, high = 10.0, size = 300)
    yp = numpy.random.uniform(low = 0.0, high = 10.0, size = 300)

    start_time = time.time()
    [indices, codes] = makeQuads(numpy.stack((xp, yp), axis = -1))

    print("Made", indices.shape[0], "quads in {0:.2f} seconds.".format(time.time() - start_time))
