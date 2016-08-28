#!/usr/bin/python
#
# Class to represent a 2D spline.
#
# Hazen 12/13
#

import math
import numpy
import numpy.linalg

import storm_analysis.sa_library.daxwriter as daxwriter

import spline1D

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
            print "Calculating spline coefficients."

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
            print "Finished calculating spline coefficients."

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


if __name__ == "__main__":

    #x = numpy.arange(0.0, 4.01, 0.1)
    #y = numpy.array([[0.0, 0.0, 0.5, 0.0, 0.0],
    #                 [0.0, 0.0, 1.5, 0.5, 0.0],
    #                 [0.0, 1.0, 2.0, 1.0, 0.0],
    #                 [0.0, 0.0, 1.5, 0.0, 0.0],
    #                 [0.0, 0.0, 0.5, 0.0, 0.0]])

    x = numpy.arange(0.0, 2.001, 0.125)
    y = numpy.array([[0.0, 1.0, 2.0],
                     [3.0, 4.0, 5.0],
                     [6.0, 7.0, 8.0]])

    s = Spline2D(y)

    surf = numpy.zeros((x.size, x.size))
    dx_surf = numpy.zeros((x.size, x.size))
    dy_surf = numpy.zeros((x.size, x.size))
    for i in range(x.size):
        for j in range(x.size):
            surf[i,j] = s.f(x[i],x[j])
            dx_surf[i,j] = s.dxf(x[i],x[j])
            dy_surf[i,j] = s.dyf(x[i],x[j])

    print surf
    #print dx_surf
    #print dy_surf

    if 0:
        surf = 100.0 + 100.0 * surf
        dx_surf = 100.0 + 100.0 * dx_surf
        dy_surf = 100.0 + 100.0 * dy_surf
        daxwriter.singleFrameDax("surf.dax", surf)
        daxwriter.singleFrameDax("dx_surf.dax", dx_surf)
        daxwriter.singleFrameDax("dy_surf.dax", dy_surf)

    #pw = pyqtgraph.plot()

    #raw_input("return to continue")    
    
