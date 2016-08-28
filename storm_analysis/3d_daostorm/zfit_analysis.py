#!/usr/bin/python
#
# (Re)do analysis from z fit onwards.
#
# Hazen 02/13
#

import sys

import library.parameters as params
import mufit_analysis

# setup
if(len(sys.argv)==3):
    mlist_file = sys.argv[1]
    parameters = params.Parameters(sys.argv[2])
else:
    print "usage: <bin> <parameters.xml>"
    exit()

alist_file = False
if(parameters.radius > 0.0):
    alist_file = mlist_file[:-9] + "alist.bin"

# z fitting
if(parameters.do_zfit):
    print "Fitting Z"
    if alist_file:
        mufit_analysis.zFitting(alist_file, parameters)
    mufit_analysis.zFitting(mlist_file, parameters)
    print ""

# drift correction
if hasattr(parameters, "drift_correction"):
    if parameters.drift_correction:
        if alist_file:
            mufit_analysis.driftCorrection([mlist_file, alist_file], parameters)
        else:
            mufit_analysis.driftCorrection([mlist_file], parameters)

#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
