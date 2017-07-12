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

def makeQuad(A,B,C,D):
    """
    Returns a quad if points A,B,C,D form a proper 
    quad, otherwise returns None.
    """

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
    dx = B[0] - A[0]
    dy = B[1] - A[1]
    dab_l = 1.0/math.sqrt(dx*dx+dy*dy)

    dx = dx * dab_l
    dy = dy * dab_l

    dab_x = [c45 * dx + s45 * dy, -s45 * dx + c45*dy]
    dab_y = [c45 * dx - s45 * dy, s45 * dx + c45*dy]

    dab_l = dab_l * root2

    # Calculate xc, yc.
    dac_x = dab_l * (C[0] - A[0])
    dac_y = dab_l * (C[1] - A[1])

    xc = dab_x[0] * dac_x + dab_x[1] * dac_y
    yc = dab_y[0] * dac_x + dab_y[1] * dac_y

    # Calcule xd, yd.
    dad_x = dab_l * (D[0] - A[0])
    dad_y = dab_l * (D[1] - A[1])

    xd = dab_x[0] * dad_x + dab_x[1] * dad_y
    yd = dab_y[0] * dad_x + dab_y[1] * dad_y

    if (xc > xd):
        return

    if ((xc + xd) > 1.0):
        return

    return MicroQuad(A,B,C,D,xc,yc,xd,yc)


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

    def isEqual(self, other, tolerance = 0.01):
        assert isinstance(other, MicroQuad)
        
        if (abs(self.xc - other.xc) > tolerance):
            return False

        if (abs(self.yc - other.yc) > tolerance):
            return False

        if (abs(self.xd - other.xd) > tolerance):
            return False

        if (abs(self.yd - other.yd) > tolerance):
            return False

        return True

        
if (__name__ == "__main__"):
    
    quad1 = makeQuad([0,0], [1,1], [0.3, 0.1], [0.6,0.1])
    quad2 = makeQuad([0,0], [2,2], [0.6, 0.2], [1.2,0.2])
    quad3 = makeQuad([0,0], [2,2], [0.5, 0.2], [1.2,0.2])
    
    print(quad1.isEqual(quad2))
    print(quad1.isEqual(quad3))
    
