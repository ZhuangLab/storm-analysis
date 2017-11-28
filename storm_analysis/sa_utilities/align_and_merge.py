#!/usr/bin/env python
"""
Align and merge two Insight3 binary files.

FIXME: This loads both files into memory so could be problematic
       if either or both are very large.

Hazen 07/17
"""

import numpy
import os
from xml.etree import ElementTree

import storm_analysis.sa_library.i3togrid as i3togrid
import storm_analysis.sa_library.imagecorrelation as imagecorrelation
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3


def alignAndMerge(file1, file2, results_file, scale = 2, dx = 0, dy = 0, z_min = -500.0, z_max = 500.0):
    assert not os.path.exists(results_file)

    z_bins = int((z_max - z_min)/50)

    # Load meta data.
    metadata1 = readinsight3.loadI3Metadata(file1)
    metadata2 = readinsight3.loadI3Metadata(file1)

    # If meta data is available, update the film length
    # field to be which ever data set is longer.
    #
    # Note that the merged file will still be messy in that the
    # frame numbers for the second movie are not changed, so they
    # will likely overlap with those of the first movie and break
    # the assumption that frame number always increases as you
    # go through the file.
    #
    if (metadata1 is not None) and (metadata2 is not None):
        f1_length = int(metadata1.find("movie").find("movie_l").text)
        f2_length = int(metadata2.find("movie").find("movie_l").text)
        if (f2_length > f1_length):
            metadata1.find("movie").find("movie_l").text = str(f2_length)
    
    i3_data1 = i3togrid.I3GData(file1, scale = scale)
    i3_data2 = i3togrid.I3GData(file2, scale = scale)

    # Determine x,y offsets.
    xy_data1 = i3_data1.i3To2DGridAllChannelsMerged()
    xy_data2 = i3_data2.i3To2DGridAllChannelsMerged()
    
    [corr, offx, offy, xy_success] = imagecorrelation.xyOffset(xy_data1,
                                                               xy_data2,
                                                               scale,
                                                               center = [dx * scale,
                                                                         dy * scale])

    assert(xy_success)

    # Update x,y positions in file2.
    offx = offx/float(scale)
    offy = offy/float(scale)
    print("x,y offsets", offx, offy)

    i3_data2.i3data['xc'] += offx
    i3_data2.i3data['yc'] += offy


    # Determine z offsets.
    if ((numpy.std(i3_data1.i3data['z']) > 1.0) and (numpy.std(i3_data2.i3data['z']) > 1.0)):
        
        xyz_data1 = i3_data1.i3To3DGridAllChannelsMerged(z_bins,
                                                         zmin = z_min,
                                                         zmax = z_max)
        xyz_data2 = i3_data2.i3To3DGridAllChannelsMerged(z_bins,
                                                         zmin = z_min,
                                                         zmax = z_max)

        [corr, fit, dz, z_success] = imagecorrelation.zOffset(xyz_data1, xyz_data2)
        assert(z_success)

        dz = dz * (z_max - z_min)/float(z_bins)
        print("z offset", dz)

        # Update z positions in file2.
        i3_data2.i3data['zc'] -= dz
    else:
        print("Data is 2D, skipping z offset calculation.")

    i3w = writeinsight3.I3Writer(results_file)
    i3w.addMolecules(i3_data1.getData())
    i3w.addMolecules(i3_data2.getData())
    if metadata1 is None:
        i3w.close()
    else:
        i3w.closeWithMetadata(ElementTree.tostring(metadata1, 'ISO-8859-1'))


if (__name__ == "__main__"):
    
    import argparse

    parser = argparse.ArgumentParser(description='Align and merge two localization lists.')

    parser.add_argument('--file1', dest='file1', type=str, required=True,
                        help = "The first localization file.")
    parser.add_argument('--file2', dest='file2', type=str, required=True,
                        help = "The second localization file.")
    parser.add_argument('--results', dest='results', type=str, required=True,
                        help = "File to save the merged results in.")
    parser.add_argument('--scale', dest='scale', type=int, required=False, default = 2,
                        help = "Scale for image cross-correlation.")
    parser.add_argument('--dx', dest='dx', type=int, required=False, default = 0,
                        help = "Initial estimate for dx.")
    parser.add_argument('--dy', dest='dy', type=int, required=False, default = 0,
                        help = "Initial estimate for dy.")
    parser.add_argument('--zmin', dest='zmin', type=float, required=False, default=-500.0,
                        help = "Minimum z value in nanometers.")
    parser.add_argument('--zmax', dest='zmax', type=float, required=False, default=500.0,
                        help = "Maximum z value in nanometers.")

    args = parser.parse_args()

    alignAndMerge(args.file1,
                  args.file2,
                  args.results,
                  scale = args.scale,
                  dx = args.dx,
                  dy = args.dy,
                  z_min = args.zmin,
                  z_max = args.zmax)

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
