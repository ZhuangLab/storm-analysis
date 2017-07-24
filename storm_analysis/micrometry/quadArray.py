#!/usr/bin/env python
"""
Uses numpy arrays instead of quad objects for greater efficiency.

Hazen 07/17
"""
import math
import numpy


def getTransform(abcd1_elt, abcd2_elt):
    """
    Returns the transform to go [A,B,C,D] in space 1 to [A,B,C,D]
    in space 2.
    """
    x = numpy.array([abcd1_elt[0], abcd1_elt[2], abcd1_elt[4], abcd1_elt[6]])
    y = numpy.array([abcd1_elt[1], abcd1_elt[3], abcd1_elt[5], abcd1_elt[7]])
    m = numpy.ones((4,3))
    for j in range(4):
        m[j,1] = abcd2_elt[2*j]
        m[j,2] = abcd2_elt[2*j+1]
        
    return [numpy.linalg.lstsq(m, x)[0],
            numpy.linalg.lstsq(m, y)[0]]
    
def makeQuadArray(kd, min_size = None, max_size = None, max_neighbors = 10):
    """
    Construct 'quads' arrays.

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

    # Access to the points in the KD tree.
    kd_data = kd.data

    # Array containing all of the possible point combinations.
    #
    # [[a1x, a1y, b1x, b1y, c1x, c1y, d1x, d1y], ..]
    #
    mn = max_neighbors
    ar_size = kd_data.shape[0]*mn*(mn-1)*(mn-1)
    ABCD = numpy.zeros((ar_size, 8))

    #
    # Iterate over points in the tree to identify groups of
    # points that could be quads.
    #
    count = 0
    for i in range(kd_data.shape[0]):
        A = kd_data[i,:]

        # Add to max_neighbors + 1 as A will always have itself as a neighbor.
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
        # Add points to ABCD array.
        #
        for j in index:
            for k in index:
                if (k == j):
                    continue
                for l in index:
                    if (l == j) or (l == k):
                        continue
                    ABCD[count,0] = A[0]
                    ABCD[count,1] = A[1]
                    ABCD[count,2] = kd_data[j,0]
                    ABCD[count,3] = kd_data[j,1]
                    ABCD[count,4] = kd_data[k,0]
                    ABCD[count,5] = kd_data[k,1]
                    ABCD[count,6] = kd_data[l,0]
                    ABCD[count,7] = kd_data[l,1]
                    count += 1

    return quadArray(ABCD[:count,:])

def quadArray(ABCD):
    """
    Construct array with quad information, as well as a reduced 
    version of the ABCD array containing only those points that 
    are actually in real quads.

    Note: Works 'in place', i.e. it modifies the input ABCD array.
    """
    c45 = math.cos(0.25 * math.pi)
    s45 = math.sin(0.25 * math.pi)
    root2 = math.sqrt(2.0)

    ar_size = ABCD.shape[0]

    # Array containing the quads.
    #
    # [[x1, y1, x2, y2], ..]
    xcyc_xdyd = numpy.zeros((ar_size, 4))

    nq = 0
    for i in range(ar_size):
        
        # Calculate scale.
        dab_x = ABCD[i,2] - ABCD[i,0]
        dab_y = ABCD[i,3] - ABCD[i,1]
        dab_l = math.sqrt(dab_x*dab_x + dab_y*dab_y)

        dab_l = 1.0/dab_l
    
        # Calculate circle center.
        cx = 0.5*(ABCD[i,0]+ABCD[i,2])
        cy = 0.5*(ABCD[i,1]+ABCD[i,3])

        # Calculate radius (squared).
        dx = ABCD[i,0] - cx
        dy = ABCD[i,1] - cy
        max_rr = dx*dx + dy*dy

        # Verify that C,D are within radius of the center point.
        dx = ABCD[i,4] - cx
        dy = ABCD[i,5] - cy
        rr = dx*dx + dy*dy
        if (rr > max_rr):
            continue

        dx = ABCD[i,6] - cx
        dy = ABCD[i,7] - cy
        rr = dx*dx + dy*dy
        if (rr > max_rr):
            continue

        # Calculate basis vectors.
        dab_x = dab_x * dab_l
        dab_y = dab_y * dab_l

        x_vec = [c45 * dab_x + s45 * dab_y, -s45 * dab_x + c45*dab_y]
        y_vec = [c45 * dab_x - s45 * dab_y, s45 * dab_x + c45*dab_y]
    
        dab_l = root2 * dab_l
        
        # Calculate xc, yc.
        dac_x = dab_l * (ABCD[i,4] - ABCD[i,0])
        dac_y = dab_l * (ABCD[i,5] - ABCD[i,1])

        xc = x_vec[0] * dac_x + x_vec[1] * dac_y
        yc = y_vec[0] * dac_x + y_vec[1] * dac_y

        # Calcule xd, yd.
        dad_x = dab_l * (ABCD[i,6] - ABCD[i,0])
        dad_y = dab_l * (ABCD[i,7] - ABCD[i,1])

        xd = x_vec[0] * dad_x + x_vec[1] * dad_y
        yd = y_vec[0] * dad_x + y_vec[1] * dad_y
    
        if (xc > xd):
            continue

        if ((xc + xd) > 1.0):
            continue

        xcyc_xdyd[nq,0] = xc
        xcyc_xdyd[nq,1] = yc
        xcyc_xdyd[nq,2] = xd
        xcyc_xdyd[nq,3] = yd

        for j in range(8):
            ABCD[nq,j] = ABCD[i,j]

        nq += 1

    return [ABCD[:nq,:], xcyc_xdyd[:nq,:]]


if (__name__ == "__main__"):

    import time

    import scipy
    import scipy.spatial

    numpy.random.seed(0)
    xp = numpy.random.uniform(low = 0.0, high = 10.0, size = 300)
    yp = numpy.random.uniform(low = 0.0, high = 10.0, size = 300)
    
    kd = scipy.spatial.KDTree(numpy.stack((xp, yp), axis = -1))

    start_time = time.time()
    [ABCD, xcyc_xdyd] = makeQuadArray(kd)

    print("Made", xcyc_xdyd.shape[0], "quads in {0:.2f} seconds.".format(time.time() - start_time))
    
    
