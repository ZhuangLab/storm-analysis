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

import sa_library.datareader as datareader
import find_peaks
import sa_library.parameters as params
import sa_library.readinsight3 as readinsight3
import sa_library.writeinsight3 as writeinsight3


src_dir = sys.path[0] + "/"


#
# Helper functions.
#

# Get "x" or "y" peak width versus z paremeters.
def getWidthParams(parameters, which, for_mu_Zfit = False):
    par = ["_wo", "_c", "_d", "A", "B", "C", "D"]
    np_par = numpy.zeros(len(par))
    for i,p in enumerate(par):
        attr = "w" + which + p
        if hasattr(parameters, attr):
            np_par[i] = getattr(parameters, attr)
    if for_mu_Zfit:
        np_par[0] = np_par[0]/parameters.pixel_size
        np_par[1] = np_par[1]*0.001
        np_par[2] = np_par[2]*0.001
#        for i in range(3):
#            np_par[i] = np_par[i]*0.001
    return np_par

# Get z range
def getZRange(parameters):
    min_z = -0.5
    max_z = 0.5
    if hasattr(parameters, "min_z"):
        min_z = parameters.min_z
    if hasattr(parameters, "max_z"):
        max_z = parameters.max_z
    return [min_z, max_z]

#
# API functions.
#

# Averages all the molecules in a track into a single molecule.
def averaging(mol_list_filename, ave_list_filename):
    proc_params = [src_dir + "../sa_utilities/avemlist",
                   mol_list_filename,
                   ave_list_filename]
    subprocess.call(proc_params)

# Performs drift correction.
def driftCorrection(list_files, parameters):
    drift_path = src_dir + "../sa_utilities/"
    drift_name = list_files[0][:-9] + "drift.txt"

    # Check if we have been asked not to do z drift correction.
    # The default is to do the correction.
    z_correct = True
    if hasattr(parameters, "z_correction"):
        if not parameters.z_correction:
            z_correct = False

    if z_correct:
        proc_params = ["python",
                       drift_path + "xyz-drift-correction.py",
                       list_files[0],
                       drift_name,
                       str(parameters.frame_step),
                       str(parameters.d_scale)]
    else:
        proc_params = ["python",
                       drift_path + "xyz-drift-correction.py",
                       list_files[0],
                       drift_name,
                       str(parameters.frame_step),
                       str(parameters.d_scale),
                       str(1)]
    subprocess.call(proc_params)

    if (os.path.exists(drift_name)):
        for list_file in list_files:
            proc_params = [drift_path + "apply-drift-correction",
                           list_file,
                           drift_name]
            subprocess.call(proc_params)

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
        wx_params = getWidthParams(parameters, "x", for_mu_Zfit = True)
        wy_params = getWidthParams(parameters, "y", for_mu_Zfit = True)
        [min_z, max_z] = getZRange(parameters)

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

# Does the frame-to-frame tracking.
def tracking(mol_list_filename, parameters):
    [min_z, max_z] = getZRange(parameters)
    proc_params = [src_dir + "../sa_utilities/tracker",
                   mol_list_filename,
                   parameters.descriptor,
                   str(parameters.radius),
                   str(1000.0*min_z),
                   str(1000.0*max_z)]
    subprocess.call(proc_params)

# Does z fitting.
def zFitting(mol_list_filename, parameters):
    if(parameters.orientation == "inverted"):
        wx_str = map(str, getWidthParams(parameters, "y"))
        wy_str = map(str, getWidthParams(parameters, "x"))
    else:
        wx_str = map(str, getWidthParams(parameters, "x"))
        wy_str = map(str, getWidthParams(parameters, "y"))
    proc_params = [src_dir + "../sa_utilities/fitz",
                   mol_list_filename, 
                   str(parameters.cutoff)] + wx_str + wy_str
    subprocess.call(proc_params)


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

    # peak finding
    print "Peak finding"
    if(not peakFinding(sys.argv[1], mlist_file, parameters)):
        print ""
        
        # tracking
        print "Tracking"
        tracking(mlist_file, parameters)

        # averaging
        alist_file = None
        if(parameters.radius > 0.0):
            alist_file = mlist_file[:-9] + "alist.bin"
            averaging(mlist_file, alist_file)
            print ""

        # z fitting
        if parameters.do_zfit:
            print "Fitting Z"
            if alist_file:
                zFitting(alist_file, parameters)
            zFitting(mlist_file, parameters)
            print ""

        # drift correction
        if hasattr(parameters, "drift_correction"):
            if parameters.drift_correction:
                print "Drift Correction"
                if alist_file:
                    driftCorrection([mlist_file, alist_file], parameters)
                else:
                    driftCorrection([mlist_file], parameters)
                print ""
    print "Analysis complete"


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
