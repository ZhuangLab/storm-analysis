#!/usr/bin/env python
"""
Measure the PSF given the raw data and localization analysis. In
theory this will be more accurate as you are using fit locations
instead of locations that were entered by hand.

2D: The expected input is a movie file that contains images
    of single molecules & the corresponding analysis file
    as created by Spliner, 3D-DAOSTORM or sCMOS.

3D: The expected input is a movie file that contains images of
    single molecules at different z positions. This movie file
    needs to have been analyzed with Spliner, 3D-DAOSTORM or sCMOS.

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
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.spliner.measure_psf_utils as measurePSFUtils


def measurePSF(movie_name, zfile_name, movie_h5_name, psf_name, want2d = False, aoi_size = 12, pixel_size = 0.1, z_range = 0.75, z_step = 0.05):
    """
    movie_name - The name of the movie file.
    zfile_name - The name of the text file containing z offset data. If this does not exist
                 then the localizations z value will be used.
    movie_h5_name - The name of the HDF5 file containing the localization information.
    psf_name - The name of the file to save the measured PSF in.
    want2d - Measure a 2D PSF.
    aoi_size - The final AOI size will 2x this number (in pixels).
    pixel_size - The pixel size in microns.
    z_range - The z range of the PSF (in microns). The actual z range is 2x z_range (i.e. 
                 from -z_range to z_range).
    z_step - The z granularity of the PSF (in microns).
    """
    # Create z scaling object.
    z_sclr = measurePSFUtils.ZScaler(z_range, z_step)
    
    # Load dax file, z offset file and molecule list file.
    dax_data = datareader.inferReader(movie_name)
    z_off = None
    if os.path.exists(zfile_name):
        data = numpy.loadtxt(zfile_name, ndmin = 2)
        valid = data[:,0]
        z_off = data[:,1]

    if want2d:
        print("Measuring 2D PSF")
    else:
        print("Measuring 3D PSF")

    # Go through the frames identifying good peaks and adding them
    # to the average psf.
    #
    max_z = z_sclr.getMaxZ()

    average_psf = numpy.zeros((max_z, 2*aoi_size, 2*aoi_size))
    peaks_used = 0
    totals = numpy.zeros(max_z, dtype = numpy.int)
    
    with saH5Py.SAH5Py(movie_h5_name) as h5:
        [dax_x, dax_y, dax_l] = dax_data.filmSize()
        for curf, locs in h5.localizationsIterator():

            # Select localizations in current frame & not near the edges.
            mask = (locs['x'] > aoi_size) & (locs['x'] < (dax_x - aoi_size - 1)) & (locs['y'] > aoi_size) & (locs['y'] < (dax_y - aoi_size - 1))
            xr = locs['y'][mask] + 1
            yr = locs['x'][mask] + 1

            # Use the z offset file if it was specified, otherwise use localization z positions.
            if z_off is None:
                if (curf == 0):
                    print("Using fit z locations.")
                zr = locs['z'][mask]
            else:
                if (curf == 0):
                    print("Using z offset file.")
                if (abs(valid[curf]) < 1.0e-6):
                    continue
                zr = numpy.ones(xr.size) * z_off[curf]

            # Remove localizations that are too close to each other.
            mask = iaUtilsC.removeNeighborsMask(xr, yr, 2.0 * aoi_size)
            print(curf, "peaks in", xr.size, ", peaks out", numpy.count_nonzero(mask))
        
            xr = xr[mask]
            yr = yr[mask]
            zr = zr[mask]

            # Use remaining localizations to calculate spline.
            image = dax_data.loadAFrame(curf).astype(numpy.float64)

            for i in range(xr.size):
                xf = xr[i]
                yf = yr[i]
                zf = zr[i]
                if want2d:
                    zi = 0
                else:
                    zi = z_sclr.convert(zf)

                # Check that the z value is in range
                if z_sclr.inRange(zi):

                    # Extract PSF.
                    psf = measurePSFUtils.extractAOI(image, aoi_size, xf, yf)

                    # Add to average psf accumulator
                    average_psf[zi,:,:] += psf
                    totals[zi] += 1

    # Check that we got at least one valid measurement.
    #
    assert (numpy.max(totals) > 0)
    
    # Set the PSF to have zero average on the X/Y boundaries.
    #
    for i in range(max_z):
        edge = numpy.concatenate((average_psf[i,0,:],
                                  average_psf[i,-1,:],
                                  average_psf[i,:,0],
                                  average_psf[i,:,-1]))
        average_psf[i,:,:] -= numpy.mean(edge)

    # Normalize the PSF.
    #
    if want2d:
        max_z = 1

    # Note: I think it makes sense to normalize to a sum of 1.0 here as the user may
    #       be using the images of single localizations as the inputs. Unlike beads
    #       we can't assume that they are all the same brightness so normalizing by
    #       the number of events would make even less sense.
    #
    for i in range(max_z):
        print("z plane {0:0d} has {1:0d} samples".format(i, totals[i]))
        if (totals[i] > 0.0):
            average_psf[i,:,:] = average_psf[i,:,:]/numpy.sum(numpy.abs(average_psf[i,:,:]))

    # Normalize to unity maximum height.
    if (numpy.max(average_psf) > 0.0):
        average_psf = average_psf/numpy.max(average_psf)
    else:
        print("Warning! Measured PSF maxima is zero or negative!")
        
    # Save PSF (in image form).
    #
    # FIXME: This may be useful but it is annoying for automated testing as this file
    #        is created in which ever directory the tests are run in.
    #
    if True:
        with tifffile.TiffWriter("psf.tif") as tf:
            for i in range(max_z):
                tf.save(average_psf[i,:,:].astype(numpy.float32))
    
    # Save PSF.
    #
    #  At least for now the PSFs use nanometers, not microns.
    #
    z_range = z_range * 1.0e+3
    z_step = z_step * 1.0e+3

    if want2d:
        psf_dict = {"psf" : average_psf[0,:,:],
                    "pixel_size" : pixel_size,
                    "type" : "2D",
                    "version" : 2.0}
        
    else:
        cur_z = -z_range
        z_vals = []
        for i in range(max_z):
            z_vals.append(cur_z)
            cur_z += z_step

        psf_dict = {"psf" : average_psf,
                    "pixel_size" : pixel_size,
                    "type" : "3D",
                    "version" : 2.0,
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
                        help = "The name of the localizations HDF5 file.")
    parser.add_argument('--psf', dest='psf', type=str, required=True,
                        help = "The name of the numpy format file to save the estimated PSF in.")
    parser.add_argument('--zoffset', dest='zoffset', type=str, required=False, default="",
                        help = "A text file with two space separated numbers on each line, the first is 1 of the frame is valid, 0 otherwise and the second is the z offset of the frame relative to the focal plane in microns.")
    parser.add_argument('--want2d', dest='want2d', action='store_true', default=False,
                        help = "Measure a 2D PSF.")
    parser.add_argument('--aoi_size', dest='aoi_size', type=int, required=False, default=12,
                        help = "The size of the area of interest around the bead in pixels. The default is 12.")
    parser.add_argument('--pixel_size', dest='pixel_size', type=float, required=False, default=100.0,
                        help = "The pixel size in nanometers. The default is 100nm.")
    parser.add_argument('--zrange', dest='zrange', type=float, required=False, default=0.750,
                        help = "The z range in microns. The PSF will be estimated from -zrange to +zrange. The default is 0.750um.")
    parser.add_argument('--zstep', dest='zstep', type=float, required=False, default=0.050,
                        help = "The z step size in microns. The default is 0.050um.")

    args = parser.parse_args()

    measurePSF(args.movie,
               args.zoffset,
               args.mlist,
               args.psf,
               want2d = args.want2d,
               aoi_size = args.aoi_size,
               pixel_size = args.pixel_size * 1.0e-3,
               z_range = args.zrange,
               z_step = args.zstep)

