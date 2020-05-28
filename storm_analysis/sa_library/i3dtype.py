#!/usr/bin/env python
"""
Insight3 data type definition & manipulation.

Hazen 4/09
"""
import numpy


def i3DataType():
    return numpy.dtype([('x', numpy.float32),   # original x location
                        ('y', numpy.float32),   # original y location
                        ('xc', numpy.float32),  # drift corrected x location
                        ('yc', numpy.float32),  # drift corrected y location
                        ('h', numpy.float32),   # fit height
                        ('a', numpy.float32),   # fit area
                        ('w', numpy.float32),   # fit width
                        ('phi', numpy.float32), # fit angle (for unconstrained elliptical gaussian)
                        ('ax', numpy.float32),  # peak aspect ratio
                        ('bg', numpy.float32),  # fit background
                        ('i', numpy.float32),   # sum - baseline for pixels included in the peak
                        ('c', numpy.int32),     # peak category ([0..9] for STORM images)
                        ('fi', numpy.int32),    # fit iterations
                        ('fr', numpy.int32),    # frame
                        ('tl', numpy.int32),    # track length
                        ('lk', numpy.int32),    # link (id of the next molecule in the trace)
                        ('z', numpy.float32),   # original z coordinate
                        ('zc', numpy.float32)]) # drift corrected z coordinate

def convertToSAHDF5(i3data, frame, nm_per_pixel):
    """
    Create a storm-analysis HDF5 compatible localizations dictionary.
    """
    peaks = convertToMultiFit(i3data, frame, nm_per_pixel)
    peaks["error"] = i3data['i']
    peaks["sum"] = i3data['a']

    return peaks

def convertToMultiFit(i3data, frame, nm_per_pixel):
    """
    Create a 3D-DAOSTORM, sCMOS or Spliner analysis compatible peak array from I3 data.

    Note that this uses the non-drift corrected positions.
    """
    i3data = maskData(i3data, (i3data['fr'] == frame))

    ax = i3data['ax']
    ww = i3data['w']
    
    peaks = {"background" : i3data['bg'],
             "height" : i3data['h'],
             "x" : i3data['x'] - 1,
             "xsigma" : 0.5*numpy.sqrt(ww*ww/ax)/nm_per_pixel,
             "y" : i3data['y'] - 1,
             "ysigma" :  0.5*numpy.sqrt(ww*ww*ax)/nm_per_pixel,
             "z" : i3data['z'] * 1.0e-3}

    return peaks

def createFromMultiFit(peaks, frame, nm_per_pixel):
    """
    Create an I3 data from the output of 3D-DAOSTORM, sCMOS or Spliner.
    """
    # Figure out how many peaks there are.
    n_peaks = peaks["x"].size

    # Create I3 data structured array.
    i3data = createDefaultI3Data(n_peaks)
    setI3Field(i3data, 'fr', frame)

    # Set fields, note X/Y axis swap.
    #
    # FIXME: Some of these properties are over-writing other properties, need
    #        to change to a new format..
    #
    if "background" in peaks:
        setI3Field(i3data, 'bg', peaks["background"])
    if "error" in peaks:
        setI3Field(i3data, 'i', peaks["error"])
    if "height" in peaks:
        setI3Field(i3data, 'h', peaks["height"])
#    if "iterations" in peaks:
#        setI3Field(i3data, 'i', peaks["iterations"])
#    if "significance" in peaks:
#        setI3Field(i3data, 'i', peaks["significance"])
    if "status" in peaks:
        setI3Field(i3data, 'fi', peaks["status"])
    if "sum" in peaks:
        setI3Field(i3data, 'a', peaks["sum"])
    if "x" in peaks:
        posSet(i3data, 'x', peaks["x"] + 1.0)
    if "xsigma" in peaks:
        wx = 2.0*peaks["xsigma"]*nm_per_pixel
        if "ysigma" in peaks:
            wy = 2.0*peaks["ysigma"]*nm_per_pixel
            setI3Field(i3data, 'ax', wy/wx)
            setI3Field(i3data, 'w', numpy.sqrt(wx*wy))
        else:
            setI3Field(i3data, 'w', wx)
    if "y" in peaks:
        posSet(i3data, 'y', peaks["y"] + 1.0)
    if "ysigma" in peaks:
        wy = 2.0*peaks["ysigma"]*nm_per_pixel
        if not ("xsigma" in peaks):
            setI3Field(i3data, 'w', wy)
    if "z" in peaks:
        posSet(i3data, 'z', peaks["z"] * 1000.0)
    if "category" in peaks:
        setI3Field(i3data, 'c', peaks["category"])

    return i3data


def createDefaultI3Data(size):
    data = numpy.zeros(size, dtype = i3DataType())
    defaults = [['x', 1.0],
                ['y', 1.0],
                ['xc', 1.0],
                ['yc', 1.0],
                ['h', 100.0],
                ['a', 10000.0],
                ['w', 300.0],
                ['phi', 0.0],
                ['ax', 1.0],
                ['bg', 0.0],
                ['i', 10000.0],
                ['c', 1],
                ['fi', 1],
                ['fr', 1],
                ['tl', 1],
                ['lk', -1],
                ['z', 0.0],
                ['zc', 0.0]]

    for elt in defaults:
        setI3Field(data, elt[0], elt[1])

    return data


def getI3DataTypeSize():
    data = numpy.zeros(1, dtype = i3DataType())
    return len(data.dtype.names)


def maskData(i3data, mask):
    """
    Creates a new i3 data structure containing only
    those elements where mask is True.
    """
    new_i3data = numpy.zeros(numpy.count_nonzero(mask), dtype = i3DataType())
    for field in i3data.dtype.names:
        new_i3data[field] = i3data[field][mask]
    return new_i3data


def posSet(i3data, field, value):
    """
    Convenience function for setting both a position
    and it's corresponding drift corrected value.
    """
    setI3Field(i3data, field, value)
    setI3Field(i3data, field + 'c', value)

    
def setI3Field(i3data, field, value):
    if field in i3data.dtype.names:
        data_type = i3data.dtype.fields[field][0]
        if isinstance(value, numpy.ndarray):
            i3data[field] = value.astype(data_type)
        else:
            i3data[field] = value * numpy.ones(i3data[field].size, dtype = data_type)


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
