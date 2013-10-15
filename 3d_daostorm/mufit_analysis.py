#!/usr/bin/python
#
# Perform mufit analysis on a dax file given parameters.
#
# Hazen 03/13
#

import numpy
import os
import random
import subprocess
import sys

import find_peaks
import sa_library.datareader as datareader
import sa_library.parameters as params
import sa_library.readinsight3 as readinsight3
import sa_library.writeinsight3 as writeinsight3
import sa_utilities.std_analysis as std_analysis

# Does the peak finding.
def peakFinding(movie_file, mlist_file, parameters):

    # open files for input & output
    movie_data = datareader.inferReader(movie_file)
    [movie_x,movie_y,movie_l] = movie_data.filmSize()

    # if the i3 file already exists, read it in,
    # write it out & start the analysis from the
    # end.
    total_peaks = 0
    if(os.path.exists(mlist_file)):
        print "Found", mlist_file
        i3data_in = readinsight3.loadI3File(mlist_file)
        try:
            curf = int(numpy.max(i3data_in['fr']))
        except ValueError:
            curf = 0
        print " Starting analysis at frame:", curf
        i3data = writeinsight3.I3Writer(mlist_file)
        if (curf > 0):
            i3data.addMolecules(i3data_in)
            total_peaks = i3data_in['x'].size
    else:
        curf = 0
        i3data = writeinsight3.I3Writer(mlist_file)

    # process parameters
    if (parameters.model == "Z"):
        wx_params = params.getWidthParams(parameters, "x", for_mu_Zfit = True)
        wy_params = params.getWidthParams(parameters, "y", for_mu_Zfit = True)
        [min_z, max_z] = params.getZRange(parameters)

        if(parameters.orientation == "inverted"):
            find_peaks.initZParams(wx_params, wy_params, min_z, max_z)
        else:
            find_peaks.initZParams(wy_params, wx_params, min_z, max_z)

    if hasattr(parameters, "start_frame"):
        if (parameters.start_frame>=curf) and (parameters.start_frame<movie_l):
            curf = parameters.start_frame

    if hasattr(parameters, "max_frame"):
        if (parameters.max_frame>0) and (parameters.max_frame<movie_l):
            movie_l = parameters.max_frame

    # analyze the movie
    # catch keyboard interrupts & "gracefully" exit.
    try:
        while(curf<movie_l):
            #for j in range(l):

            # setup analysis
            image = movie_data.loadAFrame(curf) - parameters.baseline
            mask = (image < 1.0)
            if (numpy.sum(mask) > 0):
                print " Removing negative values in frame", curf
                image[mask] = 1.0

            fdata = find_peaks.FitData(image, parameters)

            if (parameters.model == "2dfixed"):
                [peaks, residual] = find_peaks.doFit2DFixed(fdata, parameters.iterations)
            elif (parameters.model == "2d"):
                [peaks, residual] = find_peaks.doFit2D(fdata, parameters.iterations)
            elif (parameters.model == "3d"):
                [peaks, residual] = find_peaks.doFit3D(fdata, parameters.iterations)
            elif (parameters.model == "Z"):
                [peaks, residual] = find_peaks.doFitZ(fdata, parameters.iterations)
            else:
                print "Unknown model:", parameters.model

            if (type(peaks) == type(numpy.array([]))):
                # remove unconverged peaks
                peaks = find_peaks.getConvergedPeaks(peaks)

                # save results
                if(parameters.orientation == "inverted"):
                    i3data.addMultiFitMolecules(peaks, movie_x, movie_y, curf+1, parameters.pixel_size, inverted = True)
                else:
                    i3data.addMultiFitMolecules(peaks, movie_x, movie_y, curf+1, parameters.pixel_size, inverted = False)

                total_peaks += peaks.shape[0]
                print "Frame:", curf, peaks.shape[0], total_peaks
            else:
                print "Frame:", curf, 0, total_peaks
            curf += 1

        print ""
        i3data.close()
        return 0

    except KeyboardInterrupt:
        print "Analysis stopped."
        i3data.close()
        return 1


# Peform analysis if called from the command line

if __name__ == "__main__":

    # setup
    if(len(sys.argv)==3):
        parameters = params.Parameters(sys.argv[2])
        mlist_file = sys.argv[1][:-4] + "_mlist.bin"
    elif(len(sys.argv)==4):
        parameters = params.Parameters(sys.argv[3])
        mlist_file = sys.argv[2]
    else:
        print "usage: <movie> <bin> <parameters.xml>"
        exit()

    std_analysis.standardAnalysis(peakFinding,
                                  sys.argv[1],
                                  mlist_file,
                                  parameters)

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
