#!/usr/bin/env python
"""
Deconvolve images in 3D using FISTA.

This minimizes || Ax - b ||_2^2 + \lambda || x ||_1.

Hazen 11/19
"""
import numpy

import storm_analysis.sa_library.cs_decon as csDecon

import storm_analysis.fista.fista_3d as fista_3d
import storm_analysis.fista.fista_fft_c as fistaFFTC


#
# FIXME: Ignore peaks outside of user-specified AOI and/or only
#        do decon on the specified sub-region.
#
class FISTADecon(csDecon.CSDecon):

    def __init__(self, image_size, psf_object, number_zvals, timestep):
        super(FISTADecon, self).__init__(image_size, psf_object, number_zvals)

        psfs = self.createPSFs()
        
        if False:
            # Python solver (useful for debugging).
            print("Using Python solver.")
            self.cs_solver = fista_3d.FISTA(psfs, timestep)
        else:
            # C solver (about 4x faster).
            print("Using C solver.")
            self.cs_solver = fistaFFTC.FISTA(psfs, timestep)


#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
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
