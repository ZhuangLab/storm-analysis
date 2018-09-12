#!/usr/bin/env python
"""
A Python implementation of the drift algorithm described in this reference:

"Localization events-based sample drift correction for localization microscopy with redundant cross-correlation algorithm", 
Wang et al. Optics Express, 30 June 2014, Vol. 22, No. 13, DOI:10.1364/OE.22.015982.

This uses the above algorithm for XY correction, then falls back to old
approach for the Z correction.

Hazen 09/14
"""

import numpy
import pickle
import scipy.interpolate
import scipy.signal
import argparse

import storm_analysis.sa_library.drift_utilities as driftUtils
import storm_analysis.sa_library.imagecorrelation as imagecorrelation


def rccDriftCorrection(hdf5_filename, drift_filename, step, scale, z_min, z_max, correct_z, make_plots = True):
    """
    hdf5_filename - The localizations file for drift estimation.
    drift_filename - A text file to save the estimated drift in.
    step - Number of frames to group together to create a single image.
    scale - Image upsampling factor, 2.0 = 2x upsampling.
    z_min - Minimum localization z value in microns.
    z_max - Maximum localization z value in microns.
    correct_z - Estimate drift in z as well as in x/y.
    """

    z_bins = int((z_max - z_min)/0.05)
    h5_dc = driftUtils.SAH5DriftCorrection(filename = hdf5_filename,
                                           scale = scale,
                                           z_bins = z_bins)
    film_l = h5_dc.getMovieLength()
    
    max_err = 0.2

    # Sub-routines.
    def saveDriftData(fdx, fdy, fdz):
        driftUtils.saveDriftData(drift_filename, fdx, fdy, fdz)
        h5_dc.saveDriftData(fdx, fdy, fdz)

    def interpolateData(xvals, yvals):
        return driftUtils.interpolateData(xvals, yvals, film_l)

    # Don't analyze films that are empty.
    if (h5_dc.getNLocalizations() == 0):
        saveDriftData(numpy.zeros(film_l),
                      numpy.zeros(film_l),
                      numpy.zeros(film_l))
        return()
    
    # Don't analyze films that are too short.
    if (4 * step > film_l):
        saveDriftData(numpy.zeros(film_l),
                      numpy.zeros(film_l),
                      numpy.zeros(film_l))
        return

    print("Performing XY correction.")

    # Figure out how to bin the movie.
    frame = 0
    bin_edges = [0]
    while(frame < film_l):
        if ((frame + 2*step) > film_l):
            frame = film_l
        else:
            frame += step
        bin_edges.append(frame)

    # Estimate offsets between all pairs of sub images.        
    centers = []
    pairs = []
    for i in range(len(bin_edges)-1):
        centers.append((bin_edges[i+1] + bin_edges[i])/2)
        for j in range(i+1, len(bin_edges)-1):
            h5_dc.setFrameRange(bin_edges[i], bin_edges[i+1])
            sub1 = h5_dc.grid2D()

            h5_dc.setFrameRange(bin_edges[j], bin_edges[j+1])
            sub2 = h5_dc.grid2D()

            [corr, dx, dy, success] = imagecorrelation.xyOffset(sub1, sub2, scale)

            dx = dx/float(scale)
            dy = dy/float(scale)

            print("offset between frame ranges", bin_edges[i], "-", bin_edges[i+1],
                  "and", bin_edges[j], "-", bin_edges[j+1])

            if success:
                print(" -> {0:0.3f} {1:0.3f} good".format(dx, dy))
            else:
                print(" -> {0:0.3f} {1:0.3f} bad".format(dx, dy))
            print("")

            pairs.append([i, j, dx, dy, success])

    print("--")

    #
    # For testing it is faster to not have to re-run the
    # XY drift correction calculations.
    #
    #with open("test.dat", "w") as fp:
    #    pickle.dump([centers, pairs], fp)
    #
    #with open("test.dat") as fp:
    #    [centers, pairs] = pickle.load(fp)
    #

    # Prepare rij_x, rij_y, A matrix.
    rij_x = numpy.zeros(len(pairs), dtype = numpy.float32)
    rij_y = numpy.zeros(len(pairs), dtype = numpy.float32)
    A = numpy.zeros((len(pairs),len(centers)), dtype = numpy.float32)
    for i, pair in enumerate(pairs):
        rij_x[i] = pair[2]
        rij_y[i] = pair[3]
        A[i,pair[0]:pair[1]] = 1.0

    # Calculate drift (pass1). 
    # dx and dy contain the optimal offset between sub image i and sub image i+1 in x/y.
    pinv_A = numpy.linalg.pinv(A)
    dx = numpy.dot(pinv_A, rij_x)
    dy = numpy.dot(pinv_A, rij_y)

    # Calculate errors.
    err_x = numpy.dot(A, dx) - rij_x
    err_y = numpy.dot(A, dy) - rij_y

    err_d = numpy.sqrt(err_x * err_x + err_y * err_y)
    arg_sort_err = numpy.argsort(err_d)

    # Print errors before.
    if False:
        print("Before:")
        for i in range(err_d.size):
            print(i, rij_x[i], rij_y[i], A[i,:], err_d[i])
        print("")

    # Remove bad values.
    j = len(arg_sort_err) - 1

    while (j > 0) and (err_d[arg_sort_err[j]] > max_err):
        index = arg_sort_err[j]
        delA = numpy.delete(A, index, 0)
        if (numpy.linalg.matrix_rank(delA, tol = 1.0) == (len(centers)-1)):
            print(j, "removing", index, "with error", err_d[index])
            A = delA
            rij_x = numpy.delete(rij_x, index, 0)
            rij_y = numpy.delete(rij_y, index, 0)
            err_d = numpy.delete(err_d, index, 0)
            arg_sort_err[(arg_sort_err > index)] -= 1
        else:
            print("not removing", index, "with error", err_d[index])
        j -= 1

    # Print errors after.
    if False:
        print("")
        print("After:")
        for i in range(err_d.size):
            print(i, rij_x[i], rij_y[i], A[i,:], err_d[i])
        print("")

    # Calculate drift (pass2). 
    pinv_A = numpy.linalg.pinv(A)
    dx = numpy.dot(pinv_A, rij_x)
    dy = numpy.dot(pinv_A, rij_y)

    # Integrate to get final drift.
    driftx = numpy.zeros((dx.size))
    drifty = numpy.zeros((dy.size))
    for i in range(dx.size):
        driftx[i] = numpy.sum(dx[0:i])
        drifty[i] = numpy.sum(dy[0:i])

    # Print out XY results.
    for i in range(driftx.size):
        print("{0:0.1f} {1:0.3f} {2:0.3f}".format(centers[i], driftx[i], drifty[i]))

    # Create spline for interpolation.
    final_driftx = interpolateData(centers, driftx)
    final_drifty = interpolateData(centers, drifty)

    # Plot XY drift.
    if make_plots:
        import matplotlib
        import matplotlib.pyplot as pyplot

        x = numpy.arange(film_l)
        pyplot.plot(x, final_driftx, color = 'blue')
        pyplot.plot(x, final_drifty, color = 'red')
        pyplot.show()

    # Z correction.
    if not correct_z:
        saveDriftData(final_driftx,
                      final_drifty,
                      numpy.zeros(film_l))        
        h5_dc.close(verbose = False)
        return

    print("")
    print("Performing Z Correction.")

    driftz = numpy.zeros((dx.size))
    xyz_master = None
    for i in range(len(bin_edges)-1):
        h5_dc.setFrameRange(bin_edges[i], bin_edges[i+1])
        h5_dc.setDriftCorrectionXY(driftx[i], drifty[i])
        h5_dc.setDriftCorrectionZ(0.0)
        
        if xyz_master is None:
            xyz_master = h5_dc.grid3D(z_min, z_max, drift_corrected = True)
            continue

        xyz_curr = h5_dc.grid3D(z_min, z_max, drift_corrected = True)
            
        # Do z correlation
        [corr, fit, dz, z_success] = imagecorrelation.zOffset(xyz_master, xyz_curr)

        # Update Values
        if z_success:
            old_dz = dz
        else:
            dz = old_dz
            
        dz = dz * (z_max - z_min)/float(z_bins)

        if z_success:
            h5_dc.setDriftCorrectionZ(-dz)
            xyz_master += h5_dc.grid3D(z_min, z_max, drift_corrected = True)

        print("{0:d} {1:d} {2:0.3f}".format(bin_edges[i], bin_edges[i+1], dz))
        driftz[i] = -dz

    final_driftz = interpolateData(centers, driftz)

    saveDriftData(final_driftx,
                  final_drifty,
                  final_driftz)

    h5_dc.close(verbose = False)

    # Plot X,Y,Z drift.
    if make_plots:
        import matplotlib
        import matplotlib.pyplot as pyplot

        pixel_size = 160.0 # pixel size in nm.
        x = numpy.arange(film_l)
        pyplot.plot(x, pixel_size * final_driftx, color = 'red')
        pyplot.plot(x, pixel_size * final_drifty, color = 'green')
        pyplot.plot(x, 1000.0*final_driftz, color = 'blue')
        pyplot.show()


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description='Calculate drift correction following Wang, Optics Express, 2014')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "Localizations binary file to calculate drift correction from.")
    parser.add_argument('--drift', dest='drift', type=str, required=True,
                        help = "Text file to save drift correction results in.")
    parser.add_argument('--step', dest='step', type=int, required=True,
                        help = "Step size in frames.")
    parser.add_argument('--scale', dest='scale', type=int, required=True,
                        help = "Scale for up-sampled images to use for correlation. 2 is usually a good value.")
    parser.add_argument('--zmin', dest='zmin', type=float, required=False, default=-0.5,
                        help = "Minimum z value in microns.")
    parser.add_argument('--zmax', dest='zmax', type=float, required=False, default=0.5,
                        help = "Maximum z value in microns.")
    parser.add_argument('--zcorrect', dest='correct_z', type=bool, required=False, default=True,
                        help = "Also perform drift correction in Z.")

    args = parser.parse_args()

    rccDriftCorrection(args.mlist, args.drift, args.step, args.scale, args.zmin, args.zmax, args.correct_z)


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
