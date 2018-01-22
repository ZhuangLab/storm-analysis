#!/usr/bin/env python
"""
Calculate 2D FRC following Nieuwenhuizen, Nature Methods, 2013

Note that this calculates the uncorrected FRC.

Hazen 10/14
"""
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import sys

import storm_analysis
import storm_analysis.frc.frc_c as frcC
import storm_analysis.sa_library.grid_c as gridC
import storm_analysis.sa_library.sa_h5py as saH5Py


def frcCalc2d(h5_name, frc_name, scale = 8, verbose = True, show_plot = True):
    """
    Calculate the 2D FRC by putting half the tracks in the first
    image and half of them in the second image.
    """
    with saH5Py.SAH5Py(h5_name) as h5:

        pixel_size = h5.getPixelSize()
        n_tracks = h5.getNTracks()
        mx, my = h5.getMovieInformation()[:2]

        grid1 = numpy.zeros((mx * scale, my * scale), dtype = numpy.int32)
        grid2 = numpy.zeros((mx * scale, my * scale), dtype = numpy.int32)

        if verbose:
            print("Generating images.")
        for locs in h5.tracksIterator(fields = ["track_id", "x", "y"]):
            if verbose:
                sys.stdout.write(".")
                sys.stdout.flush()

            # We put all the even numbered tracks in the first image, all
            # the odd ones in the second.
            #
            mask = ((locs["track_id"]%2)==0)

            i_x = numpy.floor(locs["x"][mask] * scale).astype(numpy.int32)
            i_y = numpy.floor(locs["y"][mask] * scale).astype(numpy.int32)
            gridC.grid2D(i_x, i_y, grid1)

            i_x = numpy.floor(locs["x"][~mask] * scale).astype(numpy.int32)
            i_y = numpy.floor(locs["y"][~mask] * scale).astype(numpy.int32)
            gridC.grid2D(i_x, i_y, grid2)

        if verbose:
            sys.stdout.write("\n")

    # Compute FFT
    if verbose:
        print(numpy.sum(grid1), "points in image 1.")
        print(numpy.sum(grid2), "points in image 2.")
        print("Calculating FRC.")
        
    grid1_fft = numpy.fft.fftshift(numpy.fft.fft2(grid1))
    grid2_fft = numpy.fft.fftshift(numpy.fft.fft2(grid2))

    grid1_fft_sqr = grid1_fft * numpy.conj(grid1_fft)
    grid2_fft_sqr = grid2_fft * numpy.conj(grid2_fft)
    grid1_grid2 = grid1_fft * numpy.conj(grid2_fft)

    [frc, frc_counts] = frcC.frc(grid1_fft, grid2_fft)

    for i in range(frc.size):
        if (frc_counts[i] > 0):
            frc[i] = frc[i]/float(frc_counts[i])
        else:
            frc[i] = 0.0

    xvals = numpy.arange(frc.size)
    xvals = xvals/(float(grid1_fft.shape[0]) * pixel_size * (1.0/float(scale)))
    frc = numpy.real(frc)

    numpy.savetxt(frc_name, numpy.transpose(numpy.vstack((xvals, frc))))

    if show_plot:
        storm_analysis.configureMatplotlib()
           
        pyplot.scatter(xvals, frc, s = 4, color = 'black')
        pyplot.xlim([xvals[0], xvals[-1]])
        pyplot.ylim([-0.2,1.2])
        pyplot.xlabel("Spatial Frequency (nm-1)")
        pyplot.ylabel("Correlation")
        pyplot.show()


if (__name__ == "__main__"):

    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate 2D FRC following Nieuwenhuizen, Nature Methods, 2013')
    
    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the HDF5 localizations input file.")
    parser.add_argument('--res', dest='results', type=str, required=True,
                        help = "The name of a text file to save the results in.")
    parser.add_argument('--scale', dest='scale', type=int, required=False, default = 8,
                        help = "Scaling factor for the STORM images, default is 8.")
    parser.add_argument('--no-plot', dest='no_plot', action = 'store_true', default=False,
                        help = "Do not show a plot of the FRC curve.")

    args = parser.parse_args()

    frcCalc2d(args.hdf5, args.results, scale = args.scale, show_plot = not args.no_plot)


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
