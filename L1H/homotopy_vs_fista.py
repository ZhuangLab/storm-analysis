#!/usr/bin/python
#
# Compare performance of the homotopy versus the
# fista solver. Note that we use the lambda value
# and the l0 target from homotopy as inputs to
# the fista solver.
#
# Hazen 07/13
#

import math
import numpy
import os
import sys

import fista_lib_c
import homotopy_c
import sa_library.datareader as datareader
import sa_library.daxwriter as daxwriter
import setup_A_matrix


# Setup.

src_directory = os.path.dirname(__file__)

if(len(sys.argv)!=3):
    print "usage: homotopy_vs_fista <input_dax_file> <a_matrix>"
    exit()

in_dax_file = datareader.DaxReader(sys.argv[1])

a_mat = setup_A_matrix.loadAMatrix(sys.argv[2])

x_start = 20
y_start = 20
aoi_size = 7
keep_size = 5
scale = 8
blocks = 3

low_res = numpy.zeros((blocks*keep_size, blocks*keep_size))
ht_hres = numpy.zeros((blocks*keep_size*scale, blocks*keep_size*scale))
fista_hres = numpy.zeros((ht_hres.shape))
                     

# Initialize solvers

print "Initializing solvers"
ht = homotopy_c.Homotopy(a_mat,
                         50,
                         background_term = True,
                         positive_only = True)

fista = fista_lib_c.FISTA(a_mat,
                          positive_only = True)


# Only analyze (part of) the first frame.

ccd_baseline = 100
image = in_dax_file.loadAFrame(0) - ccd_baseline

for i in range(blocks):
    for j in range(blocks):
        lr_xo = i*keep_size + x_start
        lr_yo = j*keep_size + y_start
        lr_xf = lr_xo + aoi_size
        lr_yf = lr_yo + aoi_size
        to_analyze = image[lr_xo:lr_xf,lr_yo:lr_yf]
        eps = 1.5*math.sqrt(numpy.sum(numpy.abs(to_analyze)))

        print "Solving with homotopy", i, j
        ht.newYVector(to_analyze)
        lambda_val = ht.solve(eps, 1000)
        ht_xvec = ht.getXVector()
        l0_norm = numpy.sum(numpy.abs(ht_xvec) > 1.0)
        ht_xvec = ht_xvec[0:(scale*scale*keep_size*keep_size)]

        print "Solving with fista", lambda_val, l0_norm
        fista.newBVector(to_analyze)
        fista.iterateFISTAToL0Target(lambda_val, int(1.5*l0_norm))
        fista_xvec = fista.getXVector()[0:(scale*scale*keep_size*keep_size)]

        # save low res for reference purposes.
        lr_xo = i*keep_size
        lr_yo = j*keep_size
        lr_xf = lr_xo + keep_size
        lr_yf = lr_yo + keep_size
        low_res[lr_xo:lr_xf,lr_yo:lr_yf] = to_analyze[1:-1,1:-1]

        # save hres results.
        hr_xo = i*keep_size*scale
        hr_xf = hr_xo + keep_size*scale
        hr_yo = j*keep_size*scale
        hr_yf = hr_yo + keep_size*scale

        ht_hres[hr_xo:hr_xf,hr_yo:hr_yf] = ht_xvec.reshape((scale*keep_size,scale*keep_size))
        fista_hres[hr_xo:hr_xf,hr_yo:hr_yf] = fista_xvec.reshape((scale*keep_size,scale*keep_size))


# Save results.
ht.printProfilingData()
fista.printProfilingData()

daxwriter.singleFrameDax("low_hres.dax", low_res)
daxwriter.singleFrameDax("ht_hres.dax", ht_hres)
daxwriter.singleFrameDax("fista_hres.dax", fista_hres)


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
