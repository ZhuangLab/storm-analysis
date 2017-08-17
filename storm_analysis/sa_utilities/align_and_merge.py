#!/usr/bin/env python
"""
Align and merge two Insight3 binary files.

Hazen 07/17
"""

import numpy
import os

import storm_analysis.sa_library.i3togrid as i3togrid
import storm_analysis.sa_library.imagecorrelation as imagecorrelation
import storm_analysis.sa_library.writeinsight3 as writeinsight3


def alignAndMerge(file1, file2, results_file, scale = 2, dx = 0, dy = 0, z_min = -500.0, z_max = 500.0):
    assert not os.path.exists(results_file)

    z_bins = int((z_max - z_min)/50)
    
    i3_data1 = i3togrid.I3GData(file1, scale = scale)
    i3_data2 = i3togrid.I3GData(file2, scale = scale)

    # Determine x,y offsets.
    xy_data1 = i3_data1.i3To2DGridAllChannelsMerged()
    xy_data2 = i3_data2.i3To2DGridAllChannelsMerged()
    
    [corr, dx, dy, xy_success] = imagecorrelation.xyOffset(xy_data1,
                                                           xy_data2,
                                                           scale,
                                                           center = [dx * scale,
                                                                     dy * scale])

    assert(xy_success)

    # Update x,y positions in file2.
    dx = dx/float(scale)
    dy = dy/float(scale)
    print("x,y offsets", dx, dy)

    i3_data2.i3data['xc'] += dx
    i3_data2.i3data['yc'] += dy

    # Determine z offsets.
    xyz_data1 = i3_data1.i3To3DGridAllChannelsMerged(z_bins,
                                                     zmin = z_min,
                                                     zmax = z_max)
    xyz_data2 = i3_data2.i3To3DGridAllChannelsMerged(z_bins,
                                                     zmin = z_min,
                                                     zmax = z_max)

    [corr, fit, dz, z_success] = imagecorrelation.zOffset(xyz_data1, xyz_data2)
    assert(z_success)

    dz = dz * 1000.0/float(z_bins)
    print("z offset", dz)

    # Update z positions in file2.
    i3_data2.i3data['zc'] -= dz
    
    with writeinsight3.I3Writer(results_file) as i3w:
        i3w.addMolecules(i3_data1.getData())
        i3w.addMolecules(i3_data2.getData())
        

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
