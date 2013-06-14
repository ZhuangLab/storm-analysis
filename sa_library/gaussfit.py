#!/usr/bin/python
#
# Routines for fitting gaussians
#
# Hazen 4/09
#

import math
import numpy
import scipy
import scipy.optimize

sigma = 1.0

# Least Squares fitting
def fitAFunctionLS(data, params, fn):
    result = params
    errorfunction = lambda p: numpy.ravel(fn(*p)(*numpy.indices(data.shape)) - data)
    good = True
    [result, cov_x, infodict, mesg, success] = scipy.optimize.leastsq(errorfunction, params, full_output = 1, maxfev = 500)
    err = errorfunction(result)
    err = scipy.sum(err * err)
    if (success < 1) or (success > 4):
        print "Fitting problem!", success, mesg
        good = False
    return [result, good]

# MLE fitting, following Laurence and Chromy
def fitAFunctionMLE(data, params, fn):
    result = params
    def errorFunction(p):
        fit = fn(*p)(*numpy.indices(data.shape))
        t1 = 2.0 * numpy.sum(fit - data)
        t2 = 2.0 * numpy.sum(data * numpy.log(fit/data))
        return t1 - t2
    good = True
    try:
        [result, fopt, iter, funcalls, warnflag] = scipy.optimize.fmin(errorFunction, params, full_output = 1, maxiter = 500, disp = False)
    except:
        warnflag = 1
    if (warnflag != 0):
        print "Fitting problem!"
        good = False
    return [result, good]


#
# Fitting a gaussian, this is almost a verbatim copy of:
#
# http://www.scipy.org/Cookbook/FittingData#head-11870c917b101bb0d4b34900a0da1b7deb613bf7
#

def fixedSymmetricGaussian(background, height, center_x, center_y):
    global sigma
    width = 2.0 * sigma
    return lambda x,y: background + height*numpy.exp(-(((center_x-x)/width)**2 + ((center_y-y)/width)**2) * 2)

def symmetricGaussian1D(background, height, center_x, width):
    return lambda x: background + height*numpy.exp(-(((center_x-x)/width)**2) * 2)

def symmetricGaussian(background, height, center_x, center_y, width):
    return lambda x,y: background + height*numpy.exp(-(((center_x-x)/width)**2 + ((center_y-y)/width)**2) * 2)

def fixedEllipticalGaussian(background, height, center_x, center_y, width_x, width_y):
    return lambda x,y: background + height*numpy.exp(-(((center_x-x)/width_x)**2 + ((center_y-y)/width_y)**2) * 2)

def ellipticalGaussian(background, height, center_x, center_y, a, b, c):
    return lambda x,y: background + height*numpy.exp(-(a*(center_x-x)**2 + b*(center_x-x)*(center_y-y) + c*(center_y-y)**2))

def fitFixedSymmetricGaussian(data, a_sigma):
    # data is assumed centered on a gaussian
    # of a fixed width.
    global sigma
    sigma = a_sigma
    params = [numpy.min(data),
              numpy.max(data)-numpy.min(data),
              0.5 * data.shape[0],
              0.5 * data.shape[1]]
    return fitAFunctionLS(data, params, fixedSymmetricGaussian)

def fitFixedSymmetricGaussianMLE(data, a_sigma):
    # data is assumed centered on a gaussian
    # of a fixed width.
    global sigma
    sigma = a_sigma
    params = [numpy.min(data),
              numpy.max(data)-numpy.min(data),
              0.5 * data.shape[0],
              0.5 * data.shape[1]]
    return fitAFunctionMLE(data, params, fixedSymmetricGaussian)

def fitSymmetricGaussian1D(data, width = 0.25):
    # data is assumed centered on the gaussian
    # and of size roughly 2x the width.
    params = [numpy.min(data),
              numpy.max(data)-numpy.min(data),
              numpy.argmax(data),
              width * data.shape[0]]
    print params
    return fitAFunctionLS(data, params, symmetricGaussian1D)

def fitSymmetricGaussian(data, sigma):
    # data is assumed centered on the gaussian
    # and of size roughly 2x the width.
    params = [numpy.min(data),
              numpy.max(data),
              0.5 * data.shape[0],
              0.5 * data.shape[1],
              2.0 * sigma]
    return fitAFunctionLS(data, params, symmetricGaussian)

def fitSymmetricGaussianMLE(data, sigma):
    # data is assumed centered on the gaussian
    # and of size roughly 2x the width.
    params = [numpy.min(data),
              numpy.max(data),
              0.5 * data.shape[0],
              0.5 * data.shape[1],
              2.0 * sigma]
    return fitAFunctionMLE(data, params, symmetricGaussian)

def fitFixedEllipticalGaussian(data, sigma):
    # data is assumed centered on the gaussian
    # and of size roughly 2x the width.
    params = [numpy.min(data),
              numpy.max(data),
              0.5 * data.shape[0],
              0.5 * data.shape[1],
              2.0 * sigma,
              2.0 * sigma]
    return fitAFunctionLS(data, params, fixedEllipticalGaussian)

def fitFixedEllipticalGaussianMLE(data, sigma):
    # data is assumed centered on the gaussian
    # and of size roughly 2x the width.
    params = [numpy.min(data),
              numpy.max(data),
              0.5 * data.shape[0],
              0.5 * data.shape[1],
              2.0 * sigma,
              2.0 * sigma]
    return fitAFunctionMLE(data, params, fixedEllipticalGaussian)

def fitEllipticalGaussian(data):
    # data is assumed centered on the gaussian
    # and of size roughly 2x the width.
    params = [numpy.min(data),
              numpy.max(data),
              0.5 * data.shape[0],
              0.5 * data.shape[1],
              4.0 / data.shape[0],
              0.0,
              4.0 / data.shape[1]]
    return fitAFunctionLS(data, params, ellipticalGaussian)


#
# Fitting multiple gaussians
#
def twoSymmetricGaussian1D(height1, center_x1, width1, height2, center_x2, width2):
    return lambda x: height1*numpy.exp(-(((center_x1-x)/width1)**2) * 2) +  height2*numpy.exp(-(((center_x2-x)/width2)**2) * 2)

def fitTwoSymmetricGaussian1D(data, height1, center_x1, width1, height2, center_x2, width2):
    params = [height1,
              center_x1,
              width1,
              height2,
              center_x2,
              width2]
    return fitAFunctionLS(data, params, twoSymmetricGaussian1D)


#
# Fitting a Lorentzian
#

def lorentzian(offset, height, center, width):
    return lambda x: offset + height/3.14159 * (width/((x-center)**2 + width**2))

def fitLorentzian(data):
    params = [numpy.min(data),
              numpy.max(data) - numpy.min(data),
              numpy.argmax(data),
              1.2]
    return fitAFunctionLS(data, params, lorentzian)


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
