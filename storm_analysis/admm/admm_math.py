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

The expectation is that the contents of the cells are the FFTs of the PSFs
for different z values. This means that all the math that involves individual
PSF matrices is done element wise.

Hazen 11/19
"""
import numpy


def cellsToMatrix(cells):
    """
    Convert 2D list of matrices into a single large matrix.
    """
    nmat = len(cells)
    sx,sy = cells[0][0].shape
    m = numpy.zeros((sx*nmat, sy*nmat))
    for i in range(len(cells[0])):
        for j in range(len(cells[1])):
            m[i*sx:(i+1)*sx,j*sy:(j+1)*sy] = numpy.copy(cells[i][j])
    return m


def copyCells(cells):
    """
    Create a duplicate of a 2D list of matrices.
    """
    nmat = len(cells)
    c = emptyCells()
    for i in range(nmat):
        for j in range(nmat):
            c[i][j] = numpy.copy(cells[i][j])
    return c
    

def emptyCells(nmat):
    """
    Create an empty 2D list for storing matrices.
    """
    cells = []

    for i in range(nmat):
        row = []
        for j in range(nmat):
            row.append(None)
            
        cells.append(row)
    
    return cells


def lduG(G):
    """
    G is the 2D list of matrices containing AtA + rhoI. The A 
    matrices are the PSF matrices, I is the identity matrix and
    rho is the ADMM timestep.
    """
    nmat = len(G)
    mshape = G[0][0].shape

    # Create empty M matrix.
    M = []
    for i in range(nmat):
        row = []
        for j in range(nmat):
            row.append(numpy.zeros_like(G[0][0]))
        M.append(row)
    
    # Schur decomposition.
    D = emptyCells(nmat)
    L = emptyCells(nmat)
    U = emptyCells(nmat)

    for r in range(nmat-1,-1,-1):
        for c in range(nmat-1,-1,-1):
            k = max(r,c)
            M[c][r] = G[c][r]
            for s in range(nmat-1,k,-1):
                M[c][r] = M[c][r] - numpy.matmul(M[s][r], numpy.matmul(numpy.linalg.inv(M[s][s]), M[c][s]))
        
            if (r == c):
                D[c][r] = M[c][r]
                L[c][r] = numpy.identity(mshape[0])
                U[c][r] = numpy.identity(mshape[0])
            
            elif (r > c):
                D[c][r] = numpy.zeros(mshape)
                L[c][r] = numpy.matmul(M[c][r], numpy.linalg.inv(M[k][k]))
                U[c][r] = numpy.zeros(mshape)
            
            elif (r < c):
                D[c][r] = numpy.zeros(mshape)
                L[c][r] = numpy.zeros(mshape)
                U[c][r] = numpy.matmul(numpy.linalg.inv(M[k][k]), M[c][r])

    return [L, D, U]
            

def invD(D):
    """
    Calculate inverse of D matrices list.
    """
    nmat = len(D)
    
    d_inv = emptyCells(nmat)
    for i in range(nmat):
        for j in range(nmat):
            if (i == j):
                d_inv[i][j] = numpy.linalg.inv(D[i][j])
            else:
                d_inv[i][j] = numpy.zeros_like(D[0][0])
    return d_inv
                
                
def invL(L):
    """
    Calculate inverse of L matrices list.
    """
    nmat = len(L)
    mshape = L[0][0].shape
    
    l_tmp = copyCells(L)
    l_inv = emptyCells()
    for i in range(nmat):
        for j in range(nmat):
            if (i == j):
                l_inv[j][i] = numpy.identity(mshape[0])
            else:
                l_inv[j][i] = numpy.zeros_like(L[0][0])
                
    for j in range(nmat-1):
        for i in range(j+1,nmat):
            tmp = l_tmp[j][i]
            for k in range(nmat):
                l_tmp[k][i] = l_tmp[k][i] - numpy.matmul(tmp, l_tmp[k][j])
                l_inv[k][i] = l_inv[k][i] - numpy.matmul(tmp, l_inv[k][j])
    return l_inv
   
                
def invU(U):
    """
    Calculate inverse of U matrices list.
    """
    nmat = len(U)
    mshape = U[0][0].shape
    
    u_tmp = copyCells(U)
    u_inv = emptyCells()
    for i in range(nmat):
        for j in range(nmat):
            if (i == j):
                u_inv[j][i] = numpy.identity(mshape[0])
            else:
                u_inv[j][i] = numpy.zeros_like(U[0][0])
                
    for j in range(nmat-1,0,-1):
        for i in range(j-1,-1,-1):
            tmp = u_tmp[j][i]
            for k in range(nmat):
                u_tmp[k][i] = u_tmp[k][i] - numpy.matmul(tmp, u_tmp[k][j])
                u_inv[k][i] = u_inv[k][i] - numpy.matmul(tmp, u_inv[k][j])
    return u_inv
                

def matrixToCells(matrix, nmat):
    """
    Create list of matrices from a single large matrix.
    """
    sx = int(matrix.shape[0]/nmat)
    sy = int(matrix.shape[1]/nmat)
    
    cells = []
    for i in range(nmat):
        row = []
        for j in range(nmat):
            row.append(numpy.copy(matrix[i*sx:(i+1)*sx,j*sy:(j+1)*sy]))
        cells.append(row)
        
    return cells


def multiplyCells(A,B):
    """
    Multiply two lists of cells following matrix multiplication rules.
    """
    nmat = len(A)
    
    C = emptyCells(nmat)
    for i in range(nmat):
        for j in range(nmat):
            C[j][i] = numpy.zeros_like(A[0][0])
            for k in range(nmat):
                C[j][i] += numpy.matmul(A[j][k], B[k][i])
    return C


def printCells(cells):
    """
    Print the matrices in the cells list.
    """
    nmat = len(cells)
    
    for i in range(nmat):
        for j in range(nmat):
            print(cells[i][j])
            print()
