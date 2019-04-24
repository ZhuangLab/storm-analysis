#!/usr/bin/env python
"""
Class to represent a 1D spline.

Hazen 12/13
"""

import math
import numpy
import numpy.linalg

def roundAndCheck(x, max_x):

    if (x < 0.0):
        print("spline1D.roundAndCheck: value out of range:", x)
        return [-1, -1]
    if (x > max_x):
        print("spline1D.roundAndCheck: value out of range:", x, max_x, x - max_x)
        return [-1, -1]

    x_floor = math.floor(x)
    x_diff = x - x_floor
    ix = int(x_floor)
    if (x == max_x):
        ix -= 1
        x_diff = 1.0

    return [ix, x_diff]


class Spline(object):
    def getCoeff(self):
        return self.coeff

    def getSize(self):
        return self.max_i


class Spline1D(Spline):

    def __init__(self, y):

        self.max_i = y.size-1

        # solve for M.
        b = numpy.zeros(y.size)
        for i in range(y.size-2):
            b[i+1] = 6.0*(y[i] - 2.0*y[i+1] + y[i+2])

        A = numpy.zeros((y.size,y.size))
        A[0,0] = 1.0
        A[-1,-1] = 1.0
        for i in range(y.size-2):
            A[i+1,i] = 1.0
            A[i+1,i+1] = 4.0
            A[i+1,i+2] = 1.0

        M = numpy.linalg.solve(A,b)

        # Compute spline coefficients.
        self.coeff = numpy.zeros((y.size-1, 4))
        for i in range(y.size-1):
            self.coeff[i,3] = (M[i+1] - M[i])/6.0
            self.coeff[i,2] = M[i]/2.0
            self.coeff[i,1] = (y[i+1] - y[i]) - (M[i+1] + 2.0*M[i])/6.0
            self.coeff[i,0] = y[i]

    def dx(self, x):
        [ix, x_diff] = roundAndCheck(x, self.max_i)

        if (ix == -1):
            return 0.0

        yval = 0.0
        for i in range(3):
            yval += float(i+1) * self.coeff[ix,i+1] * math.pow(x_diff, i)
        
        return yval

    def f(self, x):
        [ix, x_diff] = roundAndCheck(x, self.max_i)

        if (ix == -1):
            return 0.0

        yval = 0.0
        for i in range(4):
            yval += self.coeff[ix,i] * math.pow(x_diff, i)
        return yval

    

if (__name__ == "__main__"):

    x = numpy.arange(0.0, 16.01, 0.1)
    y = numpy.exp(-(0.5*x-4.5)*(0.5*x-4.5))

    xv = []
    yv = []
    for i in range(x.size/10):
        xv.append(x[i*10])
        yv.append(y[i*10])
    xv.append(x[-1])
    yv.append(y[-1])
    #print xv, x.size, x[-1]
    yv = numpy.array(yv)

    s = Spline1D(yv)
    ys = numpy.zeros(x.size)
    dxs = numpy.zeros(x.size)
    for i in range(x.size):
        ys[i] = s.f(x[i])
        dxs[i] = s.dx(x[i])

    #print x[-1], dxs[0]

    import pyqtgraph
    pw = pyqtgraph.plot()
    pw.plot(x,y,pen = None, symbol = 'o')
    pw.plot(x,ys)
    #pw.plot(x,dxs)

    raw_input("return to continue")
    
