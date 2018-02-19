#!/usr/bin/env python
"""
The quad object and some helper functions.

Hazen 07/17
"""
import math
import numpy

c45 = math.cos(0.25 * math.pi)
s45 = math.sin(0.25 * math.pi)
root2 = math.sqrt(2.0)


def makeQuad(A, B, C, D):
    """
    Returns a MicroQuad if points A,B,C,D form a proper 
    quad, otherwise returns None.

    A,B define the coordinate system of the quad.
    C,D are the internal points.
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
    
    return MicroQuad(A, B, C, D, xc, yc, xd, yd)


def makeQuads(kd, min_size = None, max_size = None, max_neighbors = 10):
    """
    Construct MicroQuads.

    Note: In theory the run time of this algorithm is going to be
          proportional to the number of points times max_neighbors 
          to the 3rd power.

    kd - A scipy.spatial.KDTree object.
    min_size - A,B points must be at least this distance from each
               other.
    max_size - A,B points must be at most this distance from each 
               other.
    max_neighbors - Only consider at most this many neighbors when
               constructing quads, default is 10.
    """

    quads = []
    kd_data = kd.data

    #
    # Iterate over points in the tree to identify groups of
    # points and possibly make quads from them.
    #
    for i in range(kd_data.shape[0]):
        A = kd_data[i,:]

        # Add to max_neighbors as A will always have itself as a neighbor.
        if max_size is None:
            [dist, index] = kd.query(A, k = max_neighbors + 1)
        else:
            [dist, index] = kd.query(A, k = max_neighbors + 1, distance_upper_bound = max_size)

        #
        # Filter out points closer than the minimum distance.
        # Filter out points at infinite distance. I think these
        # are returned by KDTree when you specify both
        # 'max_neighbors' and 'distance_upper_bound'.
        #
        if min_size is None:
            mask = (dist > 1.0e-6) & (dist != numpy.inf)
        else:
            mask = (dist > min_size) & (dist != numpy.inf)

        dist = dist[mask]
        index = index[mask]

        # If we don't have at least 4 points proceed to the next A.
        if (index.size < 4):
            continue
        
        #
        # Iterate over all variations of the points in index
        # trying to create quads. Note that this is going to
        # be slow if there are lots of points to consider.
        #
        for j in index:
            B = kd_data[j,:]
            for k in index:
                if (k == j):
                    continue
                C = kd_data[k,:]
                for l in index:
                    if (l == j) or (l == k):
                        continue
                    D = kd_data[l,:]
                    quad = makeQuad(A, B, C, D)
                    if quad is not None:
                        quads.append(quad)
        
    return quads


class MicroQuad(object):

    def __init__(self, A, B, C, D, xc, yc, xd, yd):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.xc = xc
        self.xd = xd
        self.yc = yc
        self.yd = yd

    def __str__(self):
        return "Quad {0:.3f} {1:.3f} {2:.3f} {3:.3f}".format(self.xc, self.yc, self.xd, self.yd)

    def getTransform(self, other):
        """
        Returns the transform to go from self space to other space.
        """
        assert isinstance(other, MicroQuad)

        x = numpy.array([self.A[0], self.B[0], self.C[0], self.D[0]])
        y = numpy.array([self.A[1], self.B[1], self.C[1], self.D[1]])
        m = numpy.ones((4,3))
        for j, elt in enumerate([other.A, other.B, other.C, other.D]):
            m[j,1] = elt[0]
            m[j,2] = elt[1]

        return [numpy.linalg.lstsq(m, x)[0],
                numpy.linalg.lstsq(m, y)[0]]
        
    def isMatch(self, other, tolerance = 1.0e-2):
        """
        Returns True is two quads match each other.
        """
        assert isinstance(other, MicroQuad)

        #
        # There are only two ways to match:
        #
        # 1. xc1 = xc2, yc1 = yc2, xd1 = xd1, yd1 = yd2
        #
        if (abs(self.xc - other.xc) < tolerance):
            if (abs(self.yc - other.yc) < tolerance):
                if (abs(self.xd - other.xd) < tolerance):
                    if (abs(self.yd - other.yd) < tolerance):
                        return True

        #
        # 2. xc1 = yc2, yc1 = xc2, xd1 = yd1, yd1 = xd2
        #
        if (abs(self.xc - other.yc) < tolerance):
            if (abs(self.yc - other.xc) < tolerance):
                if (abs(self.xd - other.yd) < tolerance):
                    if (abs(self.yd - other.xd) < tolerance):
                        return True

        return False

        
if (__name__ == "__main__"):

    if False:
        quads = [makeQuad([0,0], [1,1], [0.3, 0.1], [0.6, 0.1]),
                 makeQuad([0,0], [1,1], [0.3, 0.1], [0.6, 0.1]),
                 makeQuad([0,0], [1,1], [0.1, 0.3], [0.1, 0.6]),
                 makeQuad([0,0], [-1,1], [-0.3, 0.1], [-0.6, 0.1]),
                 makeQuad([0,0], [1,-1], [0.3, -0.1], [0.6, -0.1]),
                 makeQuad([0,0], [-1,-1], [-0.3, -0.1], [-0.6, -0.1])]

        for i in range(1,3):
            print(quads[0].isMatch(quads[i]))
            print(quads[0].getTransform(quads[i]))
        
    if True:

        import time

        import scipy
        import scipy.spatial

        numpy.random.seed(0)
        xp = numpy.random.uniform(low = 0.0, high = 10.0, size = 300)
        yp = numpy.random.uniform(low = 0.0, high = 10.0, size = 300)
    
        kd = scipy.spatial.KDTree(numpy.stack((xp, yp), axis = -1))

        start_time = time.time()
        quads = makeQuads(kd)

        print("Made", len(quads), "quads in {0:.2f} seconds.".format(time.time() - start_time))

