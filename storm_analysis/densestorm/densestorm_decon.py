#!/usr/bin/env python
"""
Poisson sparse image deconvolution following 3DenseSTORM.

Hazen 11/19
"""
import numpy

import storm_analysis.sa_library.cs_decon as csDecon

import storm_analysis.densestorm.densestorm_3d as densestorm3D
import storm_analysis.densestorm.densestorm_c as densestormC


class DenseSTORMDecon(csDecon.CSDecon):
    
    def __init__(self, image_size, psf_object, number_zvals, beta, eta, micro):
        super(DenseSTORMDecon, self).__init__(image_size, psf_object, number_zvals)

        psfs = self.createPSFs()

        if False:
            # Python solver (useful for debugging).
            self.cs_solver = densestorm3D.DenseSTORM(psfs, beta, eta, micro)
        else:
            # C solver (faster).
            self.cs_solver = densestormC.DenseSTORM(psfs, beta, eta, micro)

    def decon(self, iterations, verbose = False):
        for i in range(iterations):
            if verbose and ((i%10) == 0):
                print(i, self.cs_solver.l2Error())
            self.cs_solver.iterate()

#
# The MIT License
#
# Copyright (c) 2019 Babcock Lab, Harvard University
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
