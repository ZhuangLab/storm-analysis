#!/usr/bin/python
#
# Simple Python interface to homotopy C library image analysis library.
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

homotopyIa = False

# C interface definition.
directory = os.path.dirname(__file__)
if not (directory == ""):
    directory += "/"

# C interface definition.
def setCInterface(homotopy_ia_lib):
    global directory
    global homotopyIa

    if(sys.platform == "win32"):
        homotopyIa = cdll.LoadLibrary(directory + homotopy_ia_lib + ".dll")
    else:
        homotopyIa = cdll.LoadLibrary(directory + homotopy_ia_lib + ".so")

    # Check that C libraries were compiled as expected.
    l1flt_size = homotopyIa.getL1FLTSize()
    if(l1flt_size != 8):
        print "L1FLT is not of type double, exiting"
        exit()

    homotopyIa.analyzeImage.argtypes = [ndpointer(dtype=numpy.float64),
                                        ndpointer(dtype=numpy.float64)]
    homotopyIa.getPeaks.argtypes = [ndpointer(dtype=numpy.float64),
                                    ndpointer(dtype=numpy.float64),
                                    ndpointer(dtype=numpy.float64),
                                    ndpointer(dtype=numpy.float64),
                                    ndpointer(dtype=numpy.int32),
                                    c_int]
    homotopyIa.getPeaks.restype = c_int
    homotopyIa.openFile.argtypes = [c_char_p]
    homotopyIa.openFile.restype = c_int
    homotopyIa.saveHighRes.argtypes = [ndpointer(dtype=numpy.float64),
                                       c_int]
    homotopyIa.saveHighRes.restype = c_int
    homotopyIa.setImageParameters.argtypes = [ndpointer(dtype=numpy.float64),
                                              c_int,
                                              c_int,
                                              c_int,
                                              c_double,
                                              c_int,
                                              c_int,
                                              c_int,
                                              c_int,
                                              c_int]

setCInterface("homotopy_ia_sse")

# Image analysis class.
class HomotopyIA:

    def __init__(self, a_matrix, epsilon, image_size, have_bg = True, positive_only = True):

        a_mat = a_matrix["a_matrix"]
        box_size = a_matrix["meas_pixels"]
        keep_size = a_matrix["keep_pixels"]
        scale = a_matrix["keep_scale"]
        overlap = (box_size - keep_size)/2

        self.image_size = image_size
        self.hr_image_size = (scale*image_size[0],scale*image_size[1])
        self.open_file = False

        print "High resolution image size:", self.hr_image_size

        # check Y size & hard exit if it is wrong.
        temp = box_size*box_size
        if(a_mat.shape[0] != temp):
            print "Unexpected number of Y elements in A matrix", a_mat.shape[0], "expected", temp
            exit()

        c_bg_term = 0
        if(have_bg):
            c_bg_term = 1

        c_a_mat = numpy.ascontiguousarray(a_mat, dtype=numpy.float64)
        
        homotopyIa.setImageParameters(c_a_mat,
                                      a_mat.shape[1],
                                      c_bg_term,
                                      positive_only,
                                      epsilon,
                                      box_size,
                                      image_size[0],
                                      image_size[1],
                                      overlap,
                                      scale)

    def analyzeImage(self, image, failure_summary = True):
        if failure_summary:
            homotopyIa.resetFailureCounter()

        c_image = numpy.ascontiguousarray(image, dtype=numpy.float64)
        c_hres_image = numpy.ascontiguousarray(numpy.zeros(self.hr_image_size), dtype=numpy.float64)
        homotopyIa.analyzeImage(c_hres_image,
                                c_image)

        if failure_summary:
            homotopyIa.printFailureCounter()

        return c_hres_image

    def closeHRDataFile(self):
        if self.open_file:
            homotopyIa.closeFile()
            self.open_file = False
        
    def findPeaks(self, image, verbose = False):
        hres_image = self.analyzeImage(image)
        return self.getPeaks(hres_image, verbose = verbose)

    def getPeaks(self, hres_image, max_peaks = 10000, verbose = False):
        c_peak_x = numpy.ascontiguousarray(numpy.zeros(max_peaks), dtype=numpy.float64)
        c_peak_y = numpy.ascontiguousarray(numpy.zeros(max_peaks), dtype=numpy.float64)
        c_peak_i = numpy.ascontiguousarray(numpy.zeros(max_peaks), dtype=numpy.float64)
        c_peak_c = numpy.ascontiguousarray(numpy.zeros(max_peaks), dtype=numpy.int32)
        num_peaks = homotopyIa.getPeaks(hres_image,
                                        c_peak_x,
                                        c_peak_y,
                                        c_peak_i,
                                        c_peak_c,
                                        max_peaks)
        if verbose:
            print "Found", num_peaks, "peaks"
        c_peak_x = c_peak_x[:num_peaks]
        c_peak_y = c_peak_y[:num_peaks]
        c_peak_i = c_peak_i[:num_peaks]
        c_peak_c = c_peak_c[:num_peaks]

        return [c_peak_x, c_peak_y, c_peak_i, c_peak_c]

    def openHRDataFile(self, file_name):
        if self.open_file:
            print "HR data file is already open."
        else:
            last_frame = homotopyIa.openFile(file_name)
            self.open_file = True
        return last_frame

    def printProfilingData(self):
        if(hasattr(homotopyIa, "printProfilingData")):
            homotopyIa.printProfilingData()

    def saveHRFrame(self, hres, frame_number):
        return homotopyIa.saveHighRes(hres, frame_number)


if (__name__ == "__main__"):

    import sys

    import sa_library.datareader as datareader
    import sa_library.parameters as parameters
    import setup_A_matrix

    if (len(sys.argv) != 3):
        print "usage <dax> <xml>"
        exit()

    dax_file = datareader.DaxReader(sys.argv[1])
    params = parameters.Parameters(sys.argv[2])

    a_mat_file = params.a_matrix

    print "Using A matrix file:", a_mat_file
    a_mat = setup_A_matrix.loadAMatrix(a_mat_file)

    image = dax_file.loadAFrame(30)
    htia = HomotopyIA(a_mat,
                      params.epsilon,
                      image.shape)

    hres_image = htia.analyzeImage(image)
    [cs_x, cs_y, cs_a, cs_i] = htia.getPeaks(hres_image)
    print "Number peaks:", cs_x.size

    htia.printProfilingData()

#
# See the accompanying license.txt file.
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
#
