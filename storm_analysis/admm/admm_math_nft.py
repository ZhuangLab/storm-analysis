#!/usr/bin/env python
"""
This module is the "no Fourier transform" version of admm_math. It exists
primarily for testing purposes.

Hazen 11/19
"""
import numpy

import storm_analysis.admm.admm_math as admmMath


def lduG(G):
    """
    G a Cell object containing AtA + rhoI matrices. The A 
    matrices are the PSF matrices, I is the identity matrix and
    rho is the ADMM timestep.
    """
    nr, nc = G.getCellsShape()
    mshape = G.getMatrixShape()

    assert (nr == nc), "G Cell must be square!"
    nmat = nr
    
    # Create empty M matrix.
    M = admmMath.Cells(nmat, nmat)
    for i in range(nmat):
        for j in range(nmat):
            M[i,j] = numpy.zeros_like(G[0,0])
    
    # Schur decomposition.
    D = admmMath.Cells(nmat, nmat)
    L = admmMath.Cells(nmat, nmat)
    U = admmMath.Cells(nmat, nmat)

    for r in range(nmat-1,-1,-1):
        for c in range(nmat-1,-1,-1):
            k = max(r,c)
            M[r,c] = G[r,c]
            for s in range(nmat-1,k,-1):
                M[r,c] = M[r,c] - numpy.matmul(M[r,s], numpy.matmul(numpy.linalg.inv(M[s,s]), M[s,c]))
        
            if (r == c):
                D[r,c] = M[r,c]
                L[r,c] = numpy.identity(mshape[0])
                U[r,c] = numpy.identity(mshape[0])
            
            elif (r > c):
                D[r,c] = numpy.zeros(mshape)
                L[r,c] = numpy.matmul(M[r,c], numpy.linalg.inv(M[k,k]))
                U[r,c] = numpy.zeros(mshape)
            
            elif (r < c):
                D[r,c] = numpy.zeros(mshape)
                L[r,c] = numpy.zeros(mshape)
                U[r,c] = numpy.matmul(numpy.linalg.inv(M[k,k]), M[r,c])

    return [L, D, U]
            

def invD(D):
    """
    Calculate inverse of D Cell.
    """
    nr, nc = D.getCellsShape()
    assert (nr == nc), "D Cell must be square!"
    nmat = nr
    
    d_inv = admmMath.Cells(nmat,nmat)
    for i in range(nmat):
        for j in range(nmat):
            if (i == j):
                d_inv[i,j] = numpy.linalg.inv(D[i,j])
            else:
                d_inv[i,j] = numpy.zeros_like(D[0,0])
    return d_inv
                
                
def invL(L):
    """
    Calculate inverse of L Cell.
    """
    nr, nc = L.getCellsShape()
    mshape = L.getMatrixShape()

    assert (nr == nc), "L Cell must be square!"
    nmat = nr
    
    l_tmp = admmMath.copyCell(L)
    l_inv = admmMath.Cells(nmat, nmat)
    for i in range(nmat):
        for j in range(nmat):
            if (i == j):
                l_inv[i,j] = numpy.identity(mshape[0])
            else:
                l_inv[i,j] = numpy.zeros_like(L[0,0])
                
    for j in range(nmat-1):
        for i in range(j+1,nmat):
            tmp = l_tmp[i,j]
            for k in range(nmat):
                l_tmp[i,k] = l_tmp[i,k] - numpy.matmul(tmp, l_tmp[j,k])
                l_inv[i,k] = l_inv[i,k] - numpy.matmul(tmp, l_inv[j,k])
    return l_inv
   
                
def invU(U):
    """
    Calculate inverse of U Cell.
    """
    nr, nc = U.getCellsShape()
    mshape = U.getMatrixShape()

    assert (nr == nc), "U Cell must be square!"
    nmat = nr
    
    u_tmp = admmMath.copyCell(U)
    u_inv = admmMath.Cells(nmat, nmat)
    for i in range(nmat):
        for j in range(nmat):
            if (i == j):
                u_inv[i,j] = numpy.identity(mshape[0])
            else:
                u_inv[i,j] = numpy.zeros_like(U[0,0])
                
    for j in range(nmat-1,0,-1):
        for i in range(j-1,-1,-1):
            tmp = u_tmp[i,j]
            for k in range(nmat):
                u_tmp[i,k] = u_tmp[i,k] - numpy.matmul(tmp, u_tmp[j,k])
                u_inv[i,k] = u_inv[i,k] - numpy.matmul(tmp, u_inv[j,k])
    return u_inv
                

def multiplyMatMat(A, B):
    """
    Multiply two Cell objects following matrix multiplication rules.
    """
    nr_a, nc_a = A.getCellsShape()
    nr_b, nc_b = B.getCellsShape()

    assert(nr_b == nc_a), "Cell sizes don't match!"
           
    C = admmMath.Cells(nr_b, nc_a)
    
    for r in range(nr_b):
        for c in range(nc_a):
            C[r,c] = numpy.zeros_like(A[0,0])
            for k in range(nr_a):
                C[r,c] += numpy.matmul(A[k,c], B[r,k])
    return C


def multiplyMatVec(A, v):
    """
    Multiply Cell object by a vector.
    """
    nr_a, nc_a = A.getCellsShape()
    mx, my = A.getMatrixShape()

    assert(v.ndim == 1), "v must be a vector!"
    assert((nr_a*my) == v.size), "A and v sizes are incorrect!"

    b = numpy.zeros((nc_a*mx))
    for r in range(nr_a):
        for c in range(nc_a):
            b[c*mx:(c+1)*mx] += numpy.matmul(A[r,c], v[r*my:(r+1)*my])

    return b


def transpose(A):
    """
    Returns the transpose of Cell.
    """
    nr, nc = A.getCellsShape()

    B = admmMath.Cells(nc, nr)
    for r in range(nr):
        for c in range(nc):
            B[c,r] = numpy.transpose(A[r,c])

    return B
