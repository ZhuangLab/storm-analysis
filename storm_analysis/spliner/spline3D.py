#!/usr/bin/python
#
# Class to represent a 3D spline.
#
# Hazen 12/13
#

import math
import numpy
import numpy.linalg

import storm_analysis.sa_library.daxwriter as daxwriter

import spline1D
import spline2D

class Spline3D(spline1D.Spline):

    def __init__(self, d, coeff = False, verbose = False):

        if (d.shape[0] != d.shape[1]) or (d.shape[0] != d.shape[2]):
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
        # Create 2D splines in the "yz-plane".
        #
        yzs = []
        for i in range(size):
            yzs.append(spline2D.Spline2D(d[i,:,:]))

        #
        # Use splines in the "yz-plane" to create splines on the
        # "x axis" with sub-integer spacing.
        #
        xs = []
        cx = 0.0
        cy = 0.0
        while (cx <= (float(self.max_i) + 0.01)):
            if (cx > float(self.max_i)):
                cx = float(self.max_i)
            cy = 0.0
            while (cy <= (float(self.max_i) + 0.01)):
                if (cy > float(self.max_i)):
                    cy = float(self.max_i)
                xv = numpy.zeros(size)
                for i in range(size):
                    xv[i] = yzs[i].f(cy,cx)
                xs.append(spline1D.Spline1D(xv))
                cy += 1.0/3.0
            cx += 1.0/3.0

        #
        # Compute spline coefficients using the "x axis" splines
        # to generate 64 values per cell and then solving for 
        # the coefficients.
        #
        self.coeff = numpy.zeros((self.max_i, self.max_i, self.max_i, 64))

        A = numpy.zeros((64,64))
        for i in range(4):
            dx = float(i)/3.0
            for j in range(4):
                dy = float(j)/3.0
                for k in range(4):
                    dz = float(k)/3.0
                    for l in range(4):
                        for m in range(4):
                            for n in range(4):
                                A[i*16+j*4+k,l*16+m*4+n] = math.pow(dx,l) * math.pow(dy,m) * math.pow(dz,n)

        b = numpy.zeros(64)
        row_size = 3*self.max_i + 1
        for i in range(self.max_i):
            for j in range(self.max_i):
                for k in range(self.max_i):
                    for m in range(4):
                        for n in range(4):
                            sp = xs[i*3*row_size + j*3 + m*row_size + n]
                            for o in range(4):
                                cx = float(k) + float(o)/3.0
                                b[m*16+n*4+o] = sp.f(cx)
                    self.coeff[i,j,k,:] = numpy.linalg.solve(A,b)

        if verbose:
            print("Finished calculating spline coefficients.")

    def dxf(self, z, y, x):
        [ix, x_diff] = spline1D.roundAndCheck(x, self.max_i)
        [iy, y_diff] = spline1D.roundAndCheck(y, self.max_i)
        [iz, z_diff] = spline1D.roundAndCheck(z, self.max_i)

        if (ix == -1) or (iy == -1) or (iz == -1):
            return 0.0

        yval = 0.0
        for i in range(3):
            for j in range(4):
                for k in range(4):
                    yval += float(i+1) * self.coeff[ix, iy, iz, (i+1)*16+j*4+k] * math.pow(x_diff, i) * math.pow(y_diff, j) * math.pow(z_diff, k)
        return yval

    def dyf(self, z, y, x):
        [ix, x_diff] = spline1D.roundAndCheck(x, self.max_i)
        [iy, y_diff] = spline1D.roundAndCheck(y, self.max_i)
        [iz, z_diff] = spline1D.roundAndCheck(z, self.max_i)

        if (ix == -1) or (iy == -1) or (iz == -1):
            return 0.0

        yval = 0.0
        for i in range(4):
            for j in range(3):
                for k in range(4):
                    yval += float(j+1) * self.coeff[ix, iy, iz, i*16+(j+1)*4+k] * math.pow(x_diff, i) * math.pow(y_diff, j) * math.pow(z_diff, k)
        return yval

    def dzf(self, z, y, x):
        [ix, x_diff] = spline1D.roundAndCheck(x, self.max_i)
        [iy, y_diff] = spline1D.roundAndCheck(y, self.max_i)
        [iz, z_diff] = spline1D.roundAndCheck(z, self.max_i)

        if (ix == -1) or (iy == -1) or (iz == -1):
            return 0.0

        yval = 0.0
        for i in range(4):
            for j in range(4):
                for k in range(3):
                    yval += float(k+1) * self.coeff[ix, iy, iz, i*16+j*4+k+1] * math.pow(x_diff, i) * math.pow(y_diff, j) * math.pow(z_diff, k)
        return yval

    def f(self, z, y, x):
        [ix, x_diff] = spline1D.roundAndCheck(x, self.max_i)
        [iy, y_diff] = spline1D.roundAndCheck(y, self.max_i)
        [iz, z_diff] = spline1D.roundAndCheck(z, self.max_i)

        if (ix == -1) or (iy == -1) or (iz == -1):
            return 0.0

        yval = 0.0
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    yval += self.coeff[ix, iy, iz, i*16+j*4+k] * math.pow(x_diff, i) * math.pow(y_diff, j) * math.pow(z_diff, k)
        return yval


if __name__ == "__main__":

    def saveStack(name, stack):
        daxf = daxwriter.DaxWriter(name, stack.shape[1], stack.shape[2])
        for i in range(stack.shape[0]):
            daxf.addFrame(stack[i,:,:])
        daxf.close()

    if 0:
        x = numpy.arange(0.0, 1.01, 1.0)
        y = numpy.array([[[0.0, 1.0],
                          [2.0, 3.0]],
                         [[4.0, 5.0],
                          [6.0, 7.0]]])
        
    if 0:
        x = numpy.arange(0.0, 2.01, 1.0)
        y = numpy.array([[[0.0, 1.0, 2.0],
                          [3.0, 4.0, 5.0],
                          [6.0, 7.0, 8.0]],
                         [[9.0, 10.0, 11.0],
                          [12.0, 13.0, 14.0],
                          [15.0, 16.0, 17.0]],
                         [[18.0, 19.0, 20.0],
                          [21.0, 22.0, 23.0],
                          [24.0, 25.0, 26.0]]])

    if 0:
        x = numpy.arange(0.0, 2.01, 0.1)
        y = numpy.array([[[0.0, 0.3, 0.0],
                          [0.0, 0.6, 0.0],
                          [0.0, 0.3, 0.0]],
                         [[0.0, 0.2, 0.0],
                          [0.2, 1.0, 0.2],
                          [0.0, 0.2, 0.0]],
                         [[0.0, 0.0, 0.0],
                          [0.3, 0.6, 0.3],
                          [0.0, 0.0, 0.0]]])

    if 1:
        x = numpy.arange(0.0, 2.01, 0.1)
        y = numpy.array([[[0.0, 1.0, 2.0],
                          [3.0, 4.0, 5.0],
                          [6.0, 7.0, 8.0]],
                         [[9.0, 10.0, 11.0],
                          [12.0, 13.0, 14.0],
                          [15.0, 16.0, 17.0]],
                         [[18.0, 19.0, 20.0],
                          [21.0, 22.0, 23.0],
                          [24.0, 25.0, 26.0]]])

    s = Spline3D(y)

    surf = numpy.zeros((x.size, x.size, x.size))
#    dx_surf = numpy.zeros((x.size, x.size, x.size))
#    dy_surf = numpy.zeros((x.size, x.size, x.size))
#    dz_surf = numpy.zeros((x.size, x.size, x.size))
    for i in range(x.size):
        for j in range(x.size):
            for k in range(x.size):
                surf[i,j,k] = s.f(x[i],x[j],x[k])
#                dx_surf[i,j,k] = s.dxf(x[i],x[j],x[k])
#                dy_surf[i,j,k] = s.dyf(x[i],x[j],x[k])
#                dz_surf[i,j,k] = s.dzf(x[i],x[j],x[k])

    if 1:
        print(surf)
#        print dx_surf
#        print dy_surf
#        print dz_surf

    if 0:
        surf = 100.0 + 100.0 * surf
        saveStack("surf.dax", surf)

    #dx_surf = 100.0 + 100.0 * dx_surf
    #dy_surf = 100.0 + 100.0 * dy_surf
    #daxwriter.singleFrameDax("surf.dax", surf)
    #daxwriter.singleFrameDax("dx_surf.dax", dx_surf)
    #daxwriter.singleFrameDax("dy_surf.dax", dy_surf)

    #pw = pyqtgraph.plot()

    #raw_input("return to continue")    
    
