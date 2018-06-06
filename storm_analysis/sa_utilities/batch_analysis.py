#!/usr/bin/env python
"""
Batch analysis for 3D-DAOSTORM, sCMOS and Spliner.

Hazen 02/14
"""

import glob
import os

import storm_analysis.sa_library.batch_run as batchRun
import storm_analysis.sa_library.datareader as datareader

def batchAnalysis(analysis_exe, input_directory, output_directory, multi_xml, max_processes = 2):
    minimum_length = 100

    # FIXME: Should also handle .tif movies?
    dax_files = glob.glob(input_directory + "*.dax")

    # Figure out which movies to analyze.
    cmd_lines = []
    for movie_file in dax_files:

        movie_obj = datareader.inferReader(movie_file)
        if(movie_obj.filmSize()[2] > minimum_length):

            print("Analyzing:", movie_file)
            basename = os.path.basename(movie_file)
            mlistname = output_directory + "/" + basename[:-4] + ".hdf5"
            cmd_lines.append(['python', analysis_exe,
                              "--movie", movie_file,
                              "--bin", mlistname,
                              "--xml", multi_xml])
    batchRun.batchRun(cmd_lines, max_processes = max_processes)

#
# The MIT License
#
# Copyright (c) 2017 Zhuang Lab, Harvard University
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


