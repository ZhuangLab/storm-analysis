#!/usr/bin/env python
"""
Merge multiple HDF5 format localization files (tracks only). No 
alignment is performed. Metadata is taken from the first file.

Hazen 01/18
"""
import sys

import storm_analysis.sa_library.sa_h5py as saH5Py


def mergeHDF5(hdf5_files, results_file):
    """
    Note: This only merges the tracks not the localizations.
    """
    with saH5Py.SAH5Py(results_file, is_existing = False) as h5_out:
        for i, h5_name in enumerate(hdf5_files):
            with saH5Py.SAH5Py(h5_name) as h5_in:
                if (i == 0):
                    [mx, my] = h5_in.getMovieInformation()[:2]
                    h5_out.setMovieInformation(mx, my, 0, "")
                    h5_out.setPixelSize(h5_in.getPixelSize())
                    h5_out.addMetadata(h5_in.getMetadata())

                for tracks in h5_in.tracksIterator():
                    sys.stdout.write(".")
                    sys.stdout.flush()
                    h5_out.addTracks(tracks)

                sys.stdout.write("\n")


if (__name__ == "__main__"):
    
    import argparse

    parser = argparse.ArgumentParser(description='Merge multiple HDF5 localization files.')

    parser.add_argument('--inbin', dest='inbin', type=str, required=True, nargs = "*",
                        help = "The names of the localization files.")
    parser.add_argument('--results', dest='results', type=str, required=True,
                        help = "File to save the merged results in.")

    args = parser.parse_args()

    mergeHDF5(args.inbin, args.results)


#
# The MIT License
#
# Copyright (c) 2018 Babcock Lab, Harvard University
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
