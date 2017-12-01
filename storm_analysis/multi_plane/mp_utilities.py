#!/usr/bin/env python
"""
Utility functions for multi-plane fitting. 

Hazen 06/17
"""
import numpy
import pickle


def getAttrs(parameters, pre, post, max_value = 8):
    pnames = []
    for i in range(max_value):
        pname = pre + str(i) + post
        if parameters.hasAttr(pname):
            pnames.append(pname)
    return pnames

def getCalibrationAttrs(parameters):
    return getAttrs(parameters, "channel", "_cal")

def getExtAttrs(parameters):
    return getAttrs(parameters, "channel", "_ext")

def getNChannels(parameters):
    return len(getExtAttrs(parameters))
    
def getOffsetAttrs(parameters):
    return getAttrs(parameters, "channel", "_offset")

def getPSFFFTAttrs(parameters):
    return getAttrs(parameters, "psf", "")

def getPupilFnAttrs(parameters):
    return getAttrs(parameters, "pupilfn", "")

def getSplineAttrs(parameters):
    return getAttrs(parameters, "spline", "")

def loadMappings(filename, margin):

    max_ch = 1
    mappings = {}

    if filename is not None:
        # Load the mappings hash.
        with open(filename, 'rb') as fp:
            mappings = pickle.load(fp)

        # Figure out number of channels.
        while "0_" + str(max_ch) + "_x" in mappings:
            max_ch += 1

    # Create and fill in mapping transform arrays.
    xt_0toN = numpy.zeros((max_ch,3))
    yt_0toN = numpy.zeros((max_ch,3))
    xt_Nto0 = numpy.zeros((max_ch,3))
    yt_Nto0 = numpy.zeros((max_ch,3))

    # 0 <-> 0 is the identity transform.
    xt_0toN[0,1] = 1.0
    yt_0toN[0,2] = 1.0
    xt_Nto0[0,1] = 1.0
    yt_Nto0[0,2] = 1.0

    for i in range(1,max_ch):
        xt_0toN[i,:] = marginCorrect(mappings["0_" + str(i) + "_x"], margin)
        yt_0toN[i,:] = marginCorrect(mappings["0_" + str(i) + "_y"], margin)
        xt_Nto0[i,:] = marginCorrect(mappings[str(i) + "_0_x"], margin)
        yt_Nto0[i,:] = marginCorrect(mappings[str(i) + "_0_y"], margin)

    xt_0toN = numpy.ascontiguousarray(xt_0toN, dtype = numpy.float64)
    yt_0toN = numpy.ascontiguousarray(yt_0toN, dtype = numpy.float64)
    xt_Nto0 = numpy.ascontiguousarray(xt_Nto0, dtype = numpy.float64)
    yt_Nto0 = numpy.ascontiguousarray(yt_Nto0, dtype = numpy.float64)
    
    return [xt_0toN, yt_0toN, xt_Nto0, yt_Nto0]
    
def marginCorrect(tr, margin):
    """
    Correct affine transform for the margin that was added to an image.
    """
    tr[0] += margin - (tr[1] + tr[2]) * margin
    return tr

