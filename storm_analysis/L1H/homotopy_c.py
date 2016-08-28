#!/usr/bin/python
#
# Simple Python interface to homotopy C library
#
# Hazen 07/12
#
# 

from ctypes import *
import math
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import sa_library.loadclib as loadclib

homotopy = False

# C interface definition.
def setCInterface(homotopy_lib):
    global homotopy

    homotopy = loadclib.loadCLibrary(os.path.dirname(__file__), homotopy_lib)

    l1flt_size = homotopy.getL1FLTSize()
    if(l1flt_size == 4):
        float_type = c_float
        float_array_type = numpy.float32
    elif(l1flt_size == 8):
        float_type = c_double
        float_array_type = numpy.float64

    homotopy.getXVector.argtypes = [ndpointer(dtype=float_array_type)]
    homotopy.initialize.argtypes = [ndpointer(dtype=float_array_type),
                                    c_int,
                                    c_int,
                                    c_int,
                                    c_int,
                                    c_int]
    homotopy.l2Error.argtypes = [ndpointer(dtype=float_array_type)]
    homotopy.l2Error.restype = float_type
    homotopy.newYVector.argtypes = [ndpointer(dtype=float_array_type)]
    homotopy.solve.argtypes = [float_type, c_int]
    homotopy.solve.restype = float_type

    if(hasattr(homotopy, "getVisited")):
        homotopy.getVisited.argtypes = [ndpointer(dtype=numpy.int32)]

    if(hasattr(homotopy, "initializeGPU")):
        homotopy.initializeGPU.argtypes = [c_char_p,
                                           c_int,
                                           c_int,
                                           c_int]

setCInterface("homotopy_storm")

# Solver class.
class Homotopy:

    def __init__(self, a_mat, max_non_zero, background_term = False, positive_only = False):

        l1flt_size = homotopy.getL1FLTSize()
        if (l1flt_size == 4):
            self.float_type = numpy.float32
        elif (l1flt_size == 8):
            self.float_type = numpy.float64

        c_a_mat = numpy.ascontiguousarray(a_mat, dtype=self.float_type)
        self.max_iters = 0
        self.nrow = c_a_mat.shape[0]
        self.ncol = c_a_mat.shape[1]
        homotopy.initialize(c_a_mat,
                            self.nrow,
                            self.ncol,
                            background_term,
                            positive_only,
                            max_non_zero)

    def cleanup(self):
        homotopy.cleanup()

    def getVisited(self):
        c_visited = numpy.ascontiguousarray(numpy.zeros(self.ncol, dtype=numpy.int32))
        homotopy.getVisited(c_visited)
        return c_visited

    def getXVector(self):
        c_x_vec = numpy.ascontiguousarray(numpy.zeros(self.ncol, dtype=self.float_type))
        homotopy.getXVector(c_x_vec)
        return c_x_vec

    def l2Error(self, x_vec):
        if (x_vector.size == self.ncol):
            c_x_vector = numpy.ascontiguousarray(x_vector, dtype=self.float_type)
            return homotopy.l2Error(c_x_vector)
        else:
            print "x vector size is not correct:", x_vector.size, "expected:", self.ncol

    def newYVector(self, y_vector):
        if (y_vector.size == self.nrow):
            c_y_vector = numpy.ascontiguousarray(y_vector, dtype=self.float_type)
            homotopy.newYVector(c_y_vector)
        else:
            print "y vector size is not correct:", y_vector.size, "expected:", self.nrow

    def printFailureCounter(self):
        homotopy.printFailureCounter()

    def printProfilingData(self):
        homotopy.printProfilingData()

    def solve(self, epsilon, max_iters):
        self.max_iters = max_iters
        self.best_lambda_val = homotopy.solve(epsilon, self.max_iters)
        return self.best_lambda_val


# GPU solver class.
class HomotopyGPU(Homotopy):

    def __init__(self, a_mat, max_non_zero, platform = 0, device = 0, work_size = 512):

        # Initialize GPU.
        fp = open(directory + "homotopy_gpu.cl")
        gpu_code = fp.read()
        homotopy.initializeGPU(gpu_code, platform, device, work_size)

        # Initialize.
        a_mat_32 = a_mat.astype(numpy.float32)
        Homotopy.__init__(self, a_mat_32, max_non_zero, False, False)


# Testing.
if __name__ == "__main__":

    A = numpy.array([[1.0,2.0,3.0],[1.0,3.0,1.5]], dtype=numpy.float64)
    y = numpy.array([6,6])

    interfaces = ["homotopy_general", "homotopy_storm", "homotopy_sse"]
    #interfaces = ["homotopy_storm", "homotopy_sse"]
    for interface in interfaces:
        print "-----", interface, "-----"
        setCInterface(interface)
        ht = Homotopy(A, 4)
        ht.newYVector(y)
        print ht.solve(0.0,4)
        print ht.getXVector()
        ht.printFailureCounter()
        ht.printProfilingData()
        print ""

#    print "---- CPU (Storm) ----"
#    setCInterface("homotopy_storm")
#    ht_storm = Homotopy(A, 4)
#    ht_storm.newYVector(y)
#    print ht_storm.solve(0.0,4)
#    print ht_storm.getXVector()
#    print ""

#    print "---- GPU ----"
#    setCInterface("homotopy_gpu")
#    ht_gpu = HomotopyGPU(A, 4)
#    ht_gpu.newYVector(y)
#    print ht_gpu.solve(0.0, 4)
#    print ht_gpu.getXVector()
#    print ""

    #print ht.getVisited()

    #
    # If this is working you should see the following:
    #
    # $ python homotopy_c.py
    # 0.0140625
    # [ 0.          1.49970703  0.99902344]
    #

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
