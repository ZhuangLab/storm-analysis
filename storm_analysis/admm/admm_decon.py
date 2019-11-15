#!/usr/bin/env python
"""
Uses ADMM to perform image deconvolution. 

This solves 1/2(Ax - b)^2 + lambda|x|.

Hazen 11/19
"""
import numpy

import storm_analysis.sa_library.cs_decon as csDecon

import storm_analysis.admm.admm_3d as admm3D
import storm_analysis.admm.admm_lasso_c as admmLassoC


class ADMMDecon(csDecon.CSDecon):
    
    def __init__(self, image_size, psf_object, number_zvals, rho):
        super(ADMMDecon, self).__init__(image_size, psf_object, number_zvals)

        psfs = self.createPSFs()

        if False:
            # Python solver (useful for debugging).
            self.cs_solver = admm3D.ADMM(psfs, rho)
        else:
            # C solver (faster).
            self.cs_solver = admmLassoC.ADMMLasso(psfs, rho)


#
# The MIT License
#
# Copyright (c) 2018 Babcock Lab, Harvard University
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
