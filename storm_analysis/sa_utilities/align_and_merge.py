#!/usr/bin/env python
"""
Align and merge two HDF5 localization files (tracks only).

Hazen 01/18
"""
import numpy
import os
import sys
import tifffile

import storm_analysis.sa_library.imagecorrelation as imagecorrelation
import storm_analysis.sa_library.sa_h5py as saH5Py


def alignAndMerge(file1, file2, results_file, scale = 2, dx = 0, dy = 0, z_min = -0.5, z_max = 0.5):
    """
    Note: This only aligns and merges the tracks not the localizations.
    """
    z_bins = int((z_max - z_min)/0.05)
    
    with saH5Py.SAH5Py(results_file, is_existing = False) as h5_out:

        # Process first file, this has no offset.
        with saH5Py.SAH5Grid(filename = file1, scale = scale, z_bins = z_bins) as h5g_in1:
            [mx, my] = h5g_in1.getMovieInformation()[:2]
            h5_out.setMovieInformation(mx, my, 0, "")
            h5_out.setPixelSize(h5g_in1.getPixelSize())
            h5_out.addMetadata(h5g_in1.getMetadata())

            for tracks in h5g_in1.tracksIterator():
                sys.stdout.write(".")
                sys.stdout.flush()
                h5_out.addTracks(tracks)

            sys.stdout.write("\n")

            im1_xy = h5g_in1.gridTracks2D()
            im1_xyz = h5g_in1.gridTracks3D(z_min, z_max)
                
        # Process second file.
        with saH5Py.SAH5Grid(filename = file2, scale = scale, z_bins = z_bins) as h5g_in2:

            # Determine X/Y offset.
            im2_xy = h5g_in2.gridTracks2D()
            [corr, offx, offy, xy_success] = imagecorrelation.xyOffset(im1_xy, im2_xy, scale,
                                                                       center = [dx * scale,
                                                                                 dy * scale])

            if False:
                tifffile.imsave("im1_xy.tif", im1_xy)
                tifffile.imsave("im2_xy.tif", im2_xy)

            assert xy_success, "Could not align images in X/Y."
            offx = offx/float(scale)
            offy = offy/float(scale)

            # Determine Z offset.
            im2_xyz = h5g_in2.gridTracks3D(z_min, z_max, dx = offx, dy = offy)

            [corr, fit, offz, z_success] = imagecorrelation.zOffset(im1_xyz, im2_xyz)
            
            assert z_success, "Could not align images in Z."
            offz = -offz * (z_max - z_min)/float(z_bins)

            for tracks in h5g_in2.tracksIterator():
                sys.stdout.write(".")
                sys.stdout.flush()
                tracks["x"] += offx
                tracks["y"] += offy
                tracks["z"] += offz
                h5_out.addTracks(tracks)

            sys.stdout.write("\n")

    return [offx, offy, offz]


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
    parser.add_argument('--zmin', dest='zmin', type=float, required=False, default=-0.5,
                        help = "Minimum z value in microns.")
    parser.add_argument('--zmax', dest='zmax', type=float, required=False, default=0.5,
                        help = "Maximum z value in microns.")

    args = parser.parse_args()

    [dx, dy, dz] = alignAndMerge(args.file1,
                                 args.file2,
                                 args.results,
                                 scale = args.scale,
                                 dx = args.dx,
                                 dy = args.dy,
                                 z_min = args.zmin,
                                 z_max = args.zmax)
    print("dx = {0:0.3f}, dy = {1:0.3f}, dz = {2:0.3f}".format(dx, dy, dz))

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
