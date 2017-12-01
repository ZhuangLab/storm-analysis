#!/usr/bin/env python
"""
Measure the PSF given the raw data and localization analysis. In
theory this will be more accurate as you are using fit locations
instead of locations that were entered by hand.

2D: The expected input is a movie file that contains images
    of single molecules & the corresponding analysis file
    as created by Spliner, 3D-DAOSTORM, sCMOS or Insight3.

3D: The expected input is a movie file that contains images of
    single molecules at different z positions. This movie file
    needs to have been analyzed with Spliner, 3D-DAOSTORM, sCMOS or 
    Insight3.

Similar to measure_psf_beads, depending on your setup you may need to change:
  1. The z range (z_range).
  2. The pixel size (pixel_size).
  3. The AOI size (aoi_size). This is less important as you
     will get to specify the final size to use when you use
     psf_to_spline.py to create the spline to use for fitting.
     It should be large enough so that there is no overlap
     between the PSFs of two different peaks.

Hazen 03/16
"""

import pickle
import numpy
import os
import scipy
import scipy.ndimage
import tifffile

import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.readinsight3 as readinsight3


def measurePSF(movie_name, zfile_name, movie_mlist, psf_name, want2d = False, aoi_size = 12, z_range = 750.0, z_step = 50.0):
    """
    The actual z range is 2x z_range (i.e. from -z_range to z_range).
    """
    
    # Load dax file, z offset file and molecule list file.
    dax_data = datareader.inferReader(movie_name)
    z_offsets = None
    if os.path.exists(zfile_name):
        try:
            z_offsets = numpy.loadtxt(zfile_name, ndmin = 2)[:,1]
        except IndexError:
            z_offsets = None
            print("z offsets were not loaded.")
    i3_data = readinsight3.loadI3File(movie_mlist)

    if want2d:
        print("Measuring 2D PSF")
    else:
        print("Measuring 3D PSF")

    #
    # Go through the frames identifying good peaks and adding them
    # to the average psf. For 3D molecule z positions are rounded to 
    # the nearest 50nm.
    #
    z_mid = int(z_range/z_step)
    max_z = 2 * z_mid + 1

    average_psf = numpy.zeros((max_z,4*aoi_size,4*aoi_size))
    peaks_used = 0
    totals = numpy.zeros(max_z)
    [dax_x, dax_y, dax_l] = dax_data.filmSize()
    for curf in range(dax_l):

        # Select localizations in current frame & not near the edges.
        mask = (i3_data['fr'] == curf+1) & (i3_data['x'] > aoi_size) & (i3_data['x'] < (dax_x - aoi_size - 1)) & (i3_data['y'] > aoi_size) & (i3_data['y'] < (dax_y - aoi_size - 1))
        xr = i3_data['x'][mask]
        yr = i3_data['y'][mask]

        # Use the z offset file if it was specified, otherwise use localization z positions.
        if z_offsets is None:
            if (curf == 0):
                print("Using fit z locations.")
            zr = i3_data['z'][mask]
        else:
            if (curf == 0):
                print("Using z offset file.")
            zr = numpy.ones(xr.size) * z_offsets[curf]

        ht = i3_data['h'][mask]

        # Remove localizations that are too close to each other.
        mask = iaUtilsC.removeNeighborsMask(xr, yr, 2.0 * aoi_size)
        print(curf, "peaks in", xr.size, ", peaks out", numpy.count_nonzero(mask))
        
        xr = xr[mask]
        yr = yr[mask]
        zr = zr[mask]
        ht = ht[mask]

        # Use remaining localizations to calculate spline.
        image = dax_data.loadAFrame(curf).astype(numpy.float64)


        for i in range(xr.size):
            xf = xr[i]
            yf = yr[i]
            zf = zr[i]
            xi = int(xf)
            yi = int(yf)
            if want2d:
                zi = 0
            else:
                zi = int(round(zf/z_step) + z_mid)

            # check the z is in range
            if (zi > -1) and (zi < max_z):

                # get localization image
                mat = image[xi-aoi_size:xi+aoi_size,
                            yi-aoi_size:yi+aoi_size]

                # zoom in by 2x
                psf = scipy.ndimage.interpolation.zoom(mat, 2.0)

                # re-center image
                psf = scipy.ndimage.interpolation.shift(psf, (-2.0*(xf-xi), -2.0*(yf-yi)), mode='nearest')

                # add to average psf accumulator
                average_psf[zi,:,:] += psf
                totals[zi] += 1

    # Force PSF to be zero (on average) at the boundaries.
    for i in range(max_z):
        edge = numpy.concatenate((average_psf[i,0,:],
                                  average_psf[i,-1,:],
                                  average_psf[i,:,0],
                                  average_psf[i,:,-1]))
        average_psf[i,:,:] -= numpy.mean(edge)

    # Normalize the PSF.
    if want2d:
        max_z = 1

    for i in range(max_z):
        print(i, totals[i])
        if (totals[i] > 0.0):
            average_psf[i,:,:] = average_psf[i,:,:]/numpy.sum(numpy.abs(average_psf[i,:,:]))

    average_psf = average_psf/numpy.max(average_psf)

    # Save PSF (in image form).
    if True:
        with tifffile.TiffWriter("psf.tif") as tf:
            for i in range(max_z):
                tf.save(average_psf[i,:,:].astype(numpy.float32))

    # Save PSF.
    if want2d:
        psf_dict = {"psf" : average_psf[0,:,:],
                    "type" : "2D"}

    else:
        cur_z = -z_range
        z_vals = []
        for i in range(max_z):
            z_vals.append(cur_z)
            cur_z += z_step

        psf_dict = {"psf" : average_psf,
                    "pixel_size" : 0.080, # 1/2 the camera pixel size in nm.
                    "type" : "3D",
                    "zmin" : -z_range,
                    "zmax" : z_range,
                    "zvals" : z_vals}

    with open(psf_name, 'wb') as fp:
        pickle.dump(psf_dict, fp)


if (__name__ == "__main__"):

    import argparse
    
    parser = argparse.ArgumentParser(description = 'Measure PSF given a movie, a list.bin file and (optionally) a z_offset file')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations file. This is a binary file in Insight3 format.")
    parser.add_argument('--psf', dest='psf', type=str, required=True,
                        help = "The name of the numpy format file to save the estimated PSF in.")
    parser.add_argument('--zoffset', dest='zoffset', type=str, required=False, default="",
                        help = "A text file with two space separated numbers on each line, the first is 1 of the frame is valid, 0 otherwise and the second is the z offset of the frame relative to the focal plane in nanometers.")
    parser.add_argument('--want2d', dest='want2d', type=bool, required=False, default=False,
                        help = "Measure a 2D PSF. The default is to measure a 3D PSF.")
    parser.add_argument('--aoi_size', dest='aoi_size', type=int, required=False, default=12,
                        help = "The size of the area of interest around the bead in pixels. The default is 12.")
    parser.add_argument('--zrange', dest='zrange', type=float, required=False, default=750,
                        help = "The z range in nanometers. The PSF will be estimated from -zrange to +zrange. The default is 750nm.")
    parser.add_argument('--zstep', dest='zstep', type=float, required=False, default=50,
                        help = "The z step size in nanometers. The default is 50nm.")

    args = parser.parse_args()

    measurePSF(args.movie, args.zoffset, args.mlist, args.psf, want2d = args.want2d, aoi_size = args.aoi_size, z_range = args.zrange, z_step = args.zstep)

