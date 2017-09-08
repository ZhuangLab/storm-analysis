#!/usr/bin/env python
"""
Merge multiple Insight3 format binary files. No alignment is 
performed. Metadata is taken from the first file.

Hazen 09/17
"""

import os
from xml.etree import ElementTree

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3


def MergeBin(bin_files, results_file):
    assert not os.path.exists(results_file)

    # Load meta data.
    metadata = readinsight3.loadI3Metadata(bin_files[0])
    
    # Create I3 writer.
    i3w = writeinsight3.I3Writer(results_file)

    # Sequentially read input files and copy into output file.
    for b_file in bin_files:
        print("Merging", b_file)
        i3_reader = readinsight3.I3Reader(b_file)
        i3_data = i3_reader.nextBlock()
        while i3_data is not False:
            print("  working..")
            i3w.addMolecules(i3_data)
            i3_data = i3_reader.nextBlock()
        
    # Close i3 output file.
    if metadata is None:
        i3w.close()
    else:
        i3w.closeWithMetadata(ElementTree.tostring(metadata, 'ISO-8859-1'))


if (__name__ == "__main__"):
    
    import argparse

    parser = argparse.ArgumentParser(description='Merge multiple localization lists.')

    parser.add_argument('--inbin', dest='inbin', type=str, required=True, nargs = "*",
                        help = "The names of the localization files.")
    parser.add_argument('--results', dest='results', type=str, required=True,
                        help = "File to save the merged results in.")

    args = parser.parse_args()

    MergeBin(args.inbin, args.results)


#
# The MIT License
#
# Copyright (c) 2017 Babcock Lab, Harvard University
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
