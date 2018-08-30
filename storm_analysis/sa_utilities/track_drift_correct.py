#!/usr/bin/env python
"""
Does tracking, z fitting and averaging and drift correction on a 
localizations binary file.

Hazen 01/18
"""
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as std_analysis


def trackDriftCorrect(h5_name, params_file):

    parameters = params.ParametersCommon().initFromFile(params_file)
    std_analysis.trackDriftCorrect(h5_name, parameters)
    

if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = '(re)does tracking, track averaging and drift correction only.')

    parser.add_argument('--bin', dest='hdf5', type=str, required=True)
    parser.add_argument('--xml', dest='settings', type=str, required=True)

    args = parser.parse_args()

    trackDriftCorrect(args.hdf5, args.settings)

#
# The MIT License
#
# Copyright (c) 2018 Zhuang Lab, Harvard University
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
