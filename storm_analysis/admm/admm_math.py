#!/usr/bin/env python
"""
This module does some of the math for doing ADMM. 

All the 3D ADMM math is a Python version of the ideas and code in the 
following references:

1. "High density 3D localization microscopy using sparse support recovery",
   Ovesny et al., Optics Express, 2014.

2. "Computational methods in single molecule localization microscopy", 
   Martin Ovesny, Thesis 2016.

3. https://github.com/zitmen/3densestorm

The expectation is that the contents of Cell objects are the FFTs of the PSFs
for different z values. This means that all the math that involves individual
PSF matrices is done element wise.

Hazen 11/19
"""
import numpy


class ADMMMathException(Exception):
    pass

##
## Cell class and cell math.
##

class Cells(object):
    """
    Class for storage and manipulation of cells of matrices 
    all of which are the same size.

    Notes:
       1. The matrices are the FFTs of PSFs, or derivatives
          thereof.

       2. The indexing is [row, col] even though internally the
          matrices are stored by column, then by row.
    """
    def __init__(self, n_rows, n_cols, mx = None, my = None, **kwds):
        """
        n_rows - Number of matrices in a row (fast axis).
        n_cols - Number of matrices in a column (slow axis).
        mx - Matrix X size (slow axis).
        my - Matrix Y size (fast axis).
        """
        super(Cells, self).__init__(**kwds)

        self.n_rows = n_rows
        self.n_cols = n_cols
        self.mx = mx
        self.my = my

        self.cells = []
        for c in range(self.n_cols):
            row = []
            for r in range(self.n_rows):
                row.append(None)
            self.cells.append(row)

    def __getitem__(self, key):
        return self.cells[key[1]][key[0]]

    def __setitem__(self, key, val):
        if (self.mx is None):
            self.mx = val.shape[0]
            if (len(val.shape) == 2):
                self.my = val.shape[1]
            
        if (self.mx != val.shape[0]):
            raise ADMMMathException("Unexpected matrix X size {0:d}, {1:d}!".format(self.mx, val.shape[0]))
        
        if (self.my is not None) and (self.my != val.shape[1]):
            raise ADMMMathException("Unexpected matrix Y size {0:d}, {1:d}!".format(self.my, val.shape[1]))

        self.cells[key[1]][key[0]] = val
        
    def getCellsShape(self):
        return (self.n_rows, self.n_cols)

    def getMatrixShape(self):
        if (self.my is None):
            return (self.mx, )
        else:
            return (self.mx, self.my)

    def getNCols(self):
        return self.n_cols

    def getNRows(self):
        return self.n_rows


def cellToMatrix(A):
    """
    Convert Cell matrices into a single large matrix.
    """
    nr, nc = A.getCellsShape()
    mx, my = A.getMatrixShape()
    
    m = numpy.zeros((nc*mx, nr*my))
    for c in range(A.getNCols()):
        for r in range(A.getNRows()):
            m[c*mx:(c+1)*mx,r*my:(r+1)*my] = numpy.copy(A[r,c])
    return m


def copyCell(A):
    """
    Create a duplicate of a Cell object.
    """
    B = Cells(*A.getCellsShape())
    for i in range(A.getNRows()):
        for j in range(A.getNCols()):
            B[i,j] = numpy.copy(A[i,j])
    return B


##
## ADMM Math
##

def lduG(G):
    """
    G a Cell object containing AtA + rhoI matrices. The A 
    matrices are the PSF matrices, I is the identity matrix and
    rho is the ADMM timestep.
    """
    nr, nc = G.getCellsShape()
    mshape = G.getMatrixShape()

    if (nr != nc):
        raise ADMMMathException("G Cell must be square!")

    nmat = nr
    
    # Create empty M matrix.
    M = Cells(nmat, nmat)
    for i in range(nmat):
        for j in range(nmat):
            M[i,j] = numpy.zeros_like(G[0,0])
    
    # Schur decomposition.
    D = Cells(nmat, nmat)
    L = Cells(nmat, nmat)
    U = Cells(nmat, nmat)

    for r in range(nmat-1,-1,-1):
        for c in range(nmat-1,-1,-1):
            k = max(r,c)
            M[r,c] = G[r,c]
            for s in range(nmat-1,k,-1):
                M[r,c] = M[r,c] - M[r,s] * M[s,c] / M[s,s]
        
            if (r == c):
                D[r,c] = M[r,c]
                L[r,c] = identityMatrix(mshape)
                U[r,c] = identityMatrix(mshape)
            
            elif (r > c):
                D[r,c] = numpy.zeros(mshape)
                L[r,c] = M[r,c] / M[k,k]
                U[r,c] = numpy.zeros(mshape)
            
            elif (r < c):
                D[r,c] = numpy.zeros(mshape)
                L[r,c] = numpy.zeros(mshape)
                U[r,c] = M[r,c] / M[k,k]

    return [L, D, U]
            

def identityMatrix(mshape, scale = 1.0):
    """
    Returns FFT of the identity matrix.
    """
    return numpy.ones(mshape, dtype = numpy.complex)*scale


def invD(D):
    """
    Calculate inverse of D Cell.
    """
    nr, nc = D.getCellsShape()

    if (nr != nc):
        raise ADMMMathException("D Cell must be square!")

    nmat = nr
    
    d_inv = Cells(nmat,nmat)
    for i in range(nmat):
        for j in range(nmat):
            if (i == j):
                d_inv[i,j] = 1.0/D[i,j]
            else:
                d_inv[i,j] = numpy.zeros_like(D[0,0])
    return d_inv
                
                
def invL(L):
    """
    Calculate inverse of L Cell.
    """
    nr, nc = L.getCellsShape()
    mshape = L.getMatrixShape()

    if (nr != nc):
        raise ADMMMathException("L Cell must be square!")

    nmat = nr
    
    l_tmp = copyCell(L)
    l_inv = Cells(nmat, nmat)
    for i in range(nmat):
        for j in range(nmat):
            if (i == j):
                l_inv[i,j] = identityMatrix(mshape)
            else:
                l_inv[i,j] = numpy.zeros_like(L[0,0])
                
    for j in range(nmat-1):
        for i in range(j+1,nmat):
            tmp = l_tmp[i,j]
            for k in range(nmat):
                l_tmp[i,k] = l_tmp[i,k] - tmp * l_tmp[j,k]
                l_inv[i,k] = l_inv[i,k] - tmp * l_inv[j,k]
    return l_inv
   
                
def invU(U):
    """
    Calculate inverse of U Cell.
    """
    nr, nc = U.getCellsShape()
    mshape = U.getMatrixShape()

    if (nr != nc):
        raise ADMMMathException("U Cell must be square!")
    
    nmat = nr
    
    u_tmp = copyCell(U)
    u_inv = Cells(nmat, nmat)
    for i in range(nmat):
        for j in range(nmat):
            if (i == j):
                u_inv[i,j] = identityMatrix(mshape)
            else:
                u_inv[i,j] = numpy.zeros_like(U[0,0])
                
    for j in range(nmat-1,0,-1):
        for i in range(j-1,-1,-1):
            tmp = u_tmp[i,j]
            for k in range(nmat):
                u_tmp[i,k] = u_tmp[i,k] - tmp * u_tmp[j,k]
                u_inv[i,k] = u_inv[i,k] - tmp * u_inv[j,k]
    return u_inv
                

def multiplyMatMat(A, B):
    """
    Multiply two Cell objects following matrix multiplication rules.
    """
    nr_a, nc_a = A.getCellsShape()
    nr_b, nc_b = B.getCellsShape()

    if (nr_b != nc_a):
        raise ADMMMathException("A, B shapes don't match!")
           
    C = Cells(nr_b, nc_a)
    
    for r in range(nr_b):
        for c in range(nc_a):
            C[r,c] = numpy.zeros_like(A[0,0])
            for k in range(nr_a):
                C[r,c] += A[k,c] * B[r,k]
    return C


def multiplyMatVec(A, v):
    """
    Multiply Cell object by a vector.
    """
    nr_a, nc_a = A.getCellsShape()
    mshape = A.getMatrixShape()

    # V is a vector.
    if (v.ndim == 1):
        mx = mshape[0]
        if (len(mshape) != 1):
            raise ADMMMathException("A and v shapes don't match!")
        if ((nr_a*mx) != v.size):
            raise ADMMMathException("A and v sizes don't match!")

        b = numpy.zeros((nc_a*mx))
        for r in range(nr_a):
            v_fft = numpy.fft.fft(v[r*mx:(r+1)*mx])
            for c in range(nc_a):
                b[c*mx:(c+1)*mx] += numpy.real(numpy.fft.ifft(A[r,c] * v_fft))

        return b

    # V is a matrix.
    #
    # In this case A must be a vector and V is the size of a single cell.
    #
    elif (v.ndim == 2):
        mx, my = mshape
        if (mx != v.shape[0]):
            raise ADMMMathException("v shape[0] doesn't match A cell size!")
        
        if (my != v.shape[1]):
            raise ADMMMathException("v shape[1] doesn't match A cell size!")

        if (nr_a != 1):
            raise ADMMMathException("A must be a vector or scalar!")
        
        b = numpy.zeros((mx, my, nc_a))
        v_fft = numpy.fft.fft2(v)
        for r in range(nr_a):
            for c in range(nc_a):
                b[:,:,c] += numpy.real(numpy.fft.ifft2(A[r,c] * v_fft))

        return b

    # V is a matrix of matrices.
    #
    # This follows the current decon convention where the first two axises
    # are the image and that last axis is the z plane, so [x, y, z].
    #
    # FIXME: When nc_a = 1, maybe we should chop the last axis?
    #
    elif (v.ndim == 3):
        mx, my = mshape
        if (mx != v.shape[0]):
            raise ADMMMathException("v shape[0] doesn't match A cell size!")
        
        if (my != v.shape[1]):
            raise ADMMMathException("v shape[1] doesn't match A cell size!")

        if (nr_a != v.shape[2]):
            raise ADMMMathException("v shape[2] doesn't match A size!")

        b = numpy.zeros((mx, my, nc_a))
        for r in range(nr_a):
            v_fft = numpy.fft.fft2(v[:,:,r])
            for c in range(nc_a):
                b[:,:,c] += numpy.real(numpy.fft.ifft2(A[r,c] * v_fft))

        return b

    else:
        raise ADMMMathException("v must be a vector or a matrix!")
                

def printCell(A):
    """
    Print the matrices in the cells list.
    """
    for r in range(A.getNRows()):
        for c in range(A.getNCols()):
            print(A[r,c])
            print()


def transpose(A):
    """
    Returns the transpose of Cell.
    """
    nr, nc = A.getCellsShape()

    B = Cells(nc, nr)
    for r in range(nr):
        for c in range(nc):
            B[c,r] = numpy.conj(A[r,c])

    return B
