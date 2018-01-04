#!/usr/bin/env python
"""
Class to represent a 2D spline.

Hazen 12/13
"""

import math
import numpy
import numpy.linalg

import storm_analysis.spliner.spline1D as spline1D

class Spline2D(spline1D.Spline):

    def __init__(self, d, coeff = False, verbose = False):

        if (d.shape[0] != d.shape[1]):
            assert "input matrix must be square!"

        size = d.shape[0]
        self.max_i = size - 1

        #
        # The coefficients have already been calculated, so just use them.
        #
        if (type(coeff) == type(numpy.array([]))):
            self.coeff = coeff
            return

        if verbose:
            print("Calculating spline coefficients.")

        #
        # Create splines along y axis.
        #
        ys = []
        for i in range(size):
            ys.append(spline1D.Spline1D(d[i,:]))

        #
        # Use splines on y axis to create splines on the
        # x axis with sub-integer spacing.
        #
        xs = []
        cx = 0.0
        while(cx <= (float(self.max_i) + 0.01)):
            if (cx > float(self.max_i)):
                cx = float(self.max_i)
            xv = numpy.zeros(size)
            for i in range(size):
                xv[i] = ys[i].f(cx)
            xs.append(spline1D.Spline1D(xv))
            cx += 1.0/3.0

        #
        # Compute spline coefficients using the x axis splines
        # to generate 16 values per grid cell and then solving
        # for the coefficients.
        #
        self.coeff = numpy.zeros((self.max_i, self.max_i, 16))

        A = numpy.zeros((16,16))
        for i in range(4):
            dx = float(i)/3.0
            for j in range(4):
                dy = float(j)/3.0
                for k in range(4):
                    for l in range(4):
                        # This is the indicing that is necessary to get d(out) to equal d(in) when printed.
                        A[i*4+j,k*4+l] = math.pow(dx,k) * math.pow(dy,l)

        b = numpy.zeros(16)
        for i in range(self.max_i):
            for j in range(self.max_i):
                for k in range(4):
                    sp = xs[3*i + k]
                    for l in range(4):
                        cy = float(j) + float(l)/3.0
                        b[k*4+l] = sp.f(cy)
                self.coeff[i,j,:] = numpy.linalg.solve(A,b)

        if verbose:
            print("Finished calculating spline coefficients.")

    def dxf(self, y, x):
        [ix, x_diff] = spline1D.roundAndCheck(x, self.max_i)
        [iy, y_diff] = spline1D.roundAndCheck(y, self.max_i)

        if (ix == -1) or (iy == -1):
            return 0.0

        yval = 0.0
        for i in range(3):
            for j in range(4):
                yval += float(i+1) * self.coeff[ix, iy, 4*(i+1)+j] * math.pow(x_diff, i) * math.pow(y_diff, j)
        return yval

    def dyf(self, y, x):
        [ix, x_diff] = spline1D.roundAndCheck(x, self.max_i)
        [iy, y_diff] = spline1D.roundAndCheck(y, self.max_i)

        if (ix == -1) or (iy == -1):
            return 0.0

        yval = 0.0
        for i in range(4):
            for j in range(3):
                yval += float(j+1) * self.coeff[ix, iy, 4*i+j+1] * math.pow(x_diff, i) * math.pow(y_diff, j)
        return yval

    def f(self, y, x):
        [ix, x_diff] = spline1D.roundAndCheck(x, self.max_i)
        [iy, y_diff] = spline1D.roundAndCheck(y, self.max_i)

        if (ix == -1) or (iy == -1):
            return 0.0

        yval = 0.0
        for i in range(4):
            for j in range(4):
                yval += self.coeff[ix, iy, 4*i+j] * math.pow(x_diff, i) * math.pow(y_diff, j)
        return yval


