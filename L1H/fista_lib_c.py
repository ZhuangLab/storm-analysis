#!/usr/bin/python
#
# Simple Python interface to fista_lib C library.
#
# Hazen 07/13
#
# 

from ctypes import *
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

# C interface definition.
directory = os.path.dirname(__file__)
if not (directory == ""):
    directory += "/"

if(sys.platform == "win32"):
    fista = cdll.LoadLibrary(directory + "fista_lib.dll")
else:
    fista = cdll.LoadLibrary(directory + "fista_lib.so")

fista.getXVector.argtypes = [ndpointer(dtype=numpy.float64)]
fista.initialize.argtypes = [ndpointer(dtype=numpy.float64),
                             c_int,
                             c_int,
                             c_int]
fista.iterateFISTA.argtypes = [c_double,
                               c_double,
                               c_int]
fista.iterateFISTAToL0Target.argtypes = [c_double,
                                         c_double,
                                         c_int]
fista.iterateFISTAToL0Target.restype = c_int
fista.iterateISTA.argtypes = [c_double,
                              c_double,
                              c_int]
fista.l1Error.restype = c_double
fista.l2Error.restype = c_double
fista.newBVector.argtypes = [ndpointer(dtype=numpy.float64)]


# Solver class.
class FISTA:

    def __init__(self, a_mat, positive_only = False):

        self.total_iters = 0

        # Calculate step size.
        norm_a = numpy.linalg.norm(a_mat)
        self.step_size = 1.0/(2.0 * norm_a * norm_a)
        
        # Initialize C library.
        c_a_mat = numpy.ascontiguousarray(a_mat, dtype = numpy.float64)
        self.nrow = c_a_mat.shape[0]
        self.ncol = c_a_mat.shape[1]
        fista.initialize(c_a_mat,
                         self.nrow,
                         self.ncol,
                         positive_only)

    def cleanup(self):
        fista.cleanup()

    def getXVector(self):
        c_x_vector = numpy.ascontiguousarray(numpy.zeros(self.ncol, dtype=numpy.float64))
        fista.getXVector(c_x_vector)
        return c_x_vector
 
    def iterateFISTA(self, lambda_term, iterations):
        self.total_iters += iterations
        fista.iterateFISTA(lambda_term,
                           self.step_size,
                           iterations)

    def iterateFISTAToL0Target(self, lambda_term, l0_target):
        iters = fista.iterateFISTAToL0Target(lambda_term,
                                             self.step_size,
                                             l0_target)
        self.total_iters += iters
        return iters

    def iterateISTA(self, lambda_term, iterations):
        self.total_iters += iterations
        fista.iterateISTA(lambda_term,
                          self.step_size,
                          iterations)

    def l1Error(self):
        return fista.l1Error()

    def l2Error(self):
        return fista.l2Error()

    def newBVector(self, b_vector):
        if (b_vector.size == self.nrow):
            c_b_vector = numpy.ascontiguousarray(b_vector, dtype=numpy.float64)
            fista.newBVector(c_b_vector)
        else:
            print "b vector size is not correct:", b_vector.size, "expected:", self.nrow

    def printProfilingData(self):
        print "Total iterations:", self.total_iters
        fista.printProfilingData()

# A quick test.
if __name__ == "__main__":

    A = numpy.array([[1.0,2.0,3.0],[1.0,3.0,1.5]], dtype=numpy.float64)
    b = numpy.array([6,6])

    solver = FISTA(A)
    solver.newBVector(b)

    for i in range(10):
        solver.iterateFISTA(0.0140625, 1)
        print i, solver.getXVector()
    
    solver.iterateFISTA(0.0140625, 300)
    print solver.getXVector()
    solver.printProfilingData()

    solver.cleanup()

#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

