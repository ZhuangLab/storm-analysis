#!/usr/bin/env python
"""
Automated XYZ drift correction for STORM movies.

Hazen 1/10

Modified to deal better with super huge insight3 files.

Hazen 11/11
"""

import numpy
import os
import scipy.signal
import sys

import storm_analysis.sa_library.arraytoimage as arraytoimag
import storm_analysis.sa_library.driftutilities as driftutilities
import storm_analysis.sa_library.i3togrid as i3togrid
import storm_analysis.sa_library.imagecorrelation as imagecorrelation


def xyzDriftCorrection(mlist_filename, drift_filename, step, scale, z_min, z_max, correct_z):
    """
    mlist_filename - The localizations file for drift estimation.
    drift_filename - A text file to save the estimated drift in.
    step - Number of frames to group together to create a single image.
    scale - Image upsampling factor, 2.0 = 2x upsampling.
    z_min - Minimum localization z value in nanometers.
    z_max - Maximum localization z value in nanoemters.
    correct_z - Estimate drift in z as well as in x/y.
    """
    
    i3_data = i3togrid.I3GDataLL(mlist_filename, scale = scale)
    film_l = i3_data.getFilmLength() - 1

    # Sub-routines.
    def saveDriftData(fdx, fdy, fdz):
        driftutilities.saveDriftData(drift_filename, fdx, fdy, fdz)

    def interpolateData(xvals, yvals):
        return driftutilities.interpolateData(xvals, yvals, film_l)

    # Don't analyze films that are too short.
    if ((4*step) >= film_l):
        saveDriftData(numpy.zeros(film_l+1),
                      numpy.zeros(film_l+1),
                      numpy.zeros(film_l+1))
        return()

    #
    # Drift correction (XY and Z are all done at the same time)
    #
    # Note that drift corrected localizations are added back into 
    # the reference image in the hopes of improving the correction
    # for subsequent localizations. 
    #

    #
    # Figure out how to bin the movie. It seemed easier to do
    # this at the beginning rather than dynamically as we
    # went through the movie.
    #
    frame = 0
    bin_edges = [0]
    while(frame < film_l):
        if ((frame + 2*step) > film_l):
            frame = film_l
        else:
            frame += step
        bin_edges.append(frame)
    
    z_bins = int((z_max - z_min)/50)
    
    xy_master = None
    xyz_master = None
    t = []
    x = []
    y = []
    z = []
    old_dx = 0.0
    old_dy = 0.0
    old_dz = 0.0
    for i in range(len(bin_edges)-1):

        # Load correct frame range.
        i3_data.loadDataInFrames(fmin = bin_edges[i], fmax = bin_edges[i+1] - 1)

        midp = (bin_edges[i+1] + bin_edges[i])/2

        xy_curr = i3_data.i3To2DGridAllChannelsMerged(uncorrected = True)

        #
        # This is to handle analysis that did not start at frame 0
        # of the movie. Basically we keep skipping ahead until we
        # find a group of frames that have some localizations.
        #
        # FIXME: There could still be problems if the movie does not
        #        start on a multiple of the step size.
        #
        if xy_master is None:
            if (numpy.sum(xy_curr) > 0):
                xy_master = xy_curr
                if correct_z:
                    xyz_master = i3_data.i3To3DGridAllChannelsMerged(z_bins,
                                                                     zmin = z_min,
                                                                     zmax = z_max,
                                                                     uncorrected = True)

            t.append(midp)
            x.append(0.0)
            y.append(0.0)
            z.append(0.0)
            print(bin_edges[i], bin_edges[i+1], numpy.sum(xy_curr), 0.0, 0.0, 0.0)
            continue
                
        # Correlate to master image.
        [corr, dx, dy, xy_success] = imagecorrelation.xyOffset(xy_master,
                                                               xy_curr,
                                                               i3_data.getScale(),
                                                               center = [x[i-1] * scale,
                                                                         y[i-1] * scale])

        #
        # Update values. If we failed, we just use the last successful
        # offset measurement and hope this is close enough.
        #
        if xy_success:
            old_dx = dx
            old_dy = dy
        else:
            dx = old_dx
            dy = old_dy

        dx = dx/float(scale)
        dy = dy/float(scale)

        t.append(midp)
        x.append(dx)
        y.append(dy)

        #
        # Apply the x/y drift correction to the current 'test'
        # localizations and add them into the master, but only
        # if the offset was measured successfully.
        #
        i3_data.applyXYDriftCorrection(dx,dy)
        if xy_success:
            # Add current to master
            xy_master += i3_data.i3To2DGridAllChannelsMerged()

        #
        # Do Z correlation if requested.
        #
        dz = old_dz
        if correct_z and xy_success:

            xyz_curr = i3_data.i3To3DGridAllChannelsMerged(z_bins,
                                                           zmin = z_min,
                                                           zmax = z_max,
                                                           uncorrected = True)

            # Do z correlation
            [corr, fit, dz, z_success] = imagecorrelation.zOffset(xyz_master, xyz_curr)

            # Update Values
            if z_success:
                old_dz = dz
            else:
                dz = old_dz
            
            dz = dz * (z_max - z_min)/float(z_bins)

            if z_success:
                i3_data.applyZDriftCorrection(-dz)
                xyz_master += i3_data.i3To3DGridAllChannelsMerged(z_bins,
                                                                  zmin = z_min,
                                                                  zmax = z_max)

        z.append(dz)

        print(bin_edges[i], bin_edges[i+1], numpy.sum(xy_curr), dx, dy, dz)

    i3_data.close()

    #
    # Create numpy versions of the drift arrays. We estimated the drift
    # for groups of frames. We use interpolation to create an estimation
    # for each individual frame.
    #
    nt = numpy.array(t)
    final_driftx = interpolateData(nt, numpy.array(x))
    final_drifty = interpolateData(nt, numpy.array(y))
    final_driftz = interpolateData(nt, numpy.array(z))

    saveDriftData(final_driftx,
                  final_drifty,
                  final_driftz)

    
if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description='Calculate drift correction using image correlation')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "Localizations binary file to calculate drift correction from.")
    parser.add_argument('--drift', dest='drift', type=str, required=True,
                        help = "Text file to save drift correction results in.")
    parser.add_argument('--step', dest='step', type=int, required=True,
                        help = "Step size in frames.")
    parser.add_argument('--scale', dest='scale', type=int, required=True,
                        help = "Scale for up-sampled images to use for correlation. 2 is usually a good value.")
    parser.add_argument('--zmin', dest='zmin', type=float, required=False, default=-500.0,
                        help = "Minimum z value in nanometers.")
    parser.add_argument('--zmax', dest='zmax', type=float, required=False, default=500.0,
                        help = "Maximum z value in nanometers.")
    parser.add_argument('--zcorrect', dest='correct_z', type=bool, required=False, default=True,
                        help = "Also perform drift correction in Z.")

    args = parser.parse_args()

    xyzDriftCorrection(args.mlist, args.drift, args.step, args.scale, args.zmin, args.zmax, args.correct_z)

    
#
# The MIT License
#
# Copyright (c) 2014 Zhuang Lab, Harvard University
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
