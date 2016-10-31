#!/usr/bin/python
#
# Measure the 3D PSF a movie given the locations of the beads
# of interest in the movie and the z-offset of each frame of
# the movie. It is assumed that the drift over the time
# course of the movie is neglible.
#
# Depending on your setup you may need to change:
#  1. The z range (z_range).
#  2. The pixel size (pixel_size).
#  3. The AOI size (aoi_size). This is less important as you
#     will get to specify the final size to use when you use
#     psf_to_spline.py to create the spline to use for fitting.
#
# Hazen 1/16
#

import pickle
import numpy
import scipy
import scipy.ndimage
import sys

import storm_analysis.sa_library.datareader as datareader

def measurePSFBeads(movie_name, zfile_name, beads_file, psf_name, want2d = False, aoi_size = 12, z_range = 600.0, z_step = 50.0):

    # Load movie file.
    movie_data = datareader.inferReader(movie_name)

    #
    # Load the z-offset information for the dax file.
    #
    #   This is a text file with one line per frame that contains the 
    #   z-offset (in nm) for that frame. Each line is a space separated
    #   valid, z_pos pair. If valid if 0 the frame will be ignored,
    #   otherwise it will be used.
    #
    data = numpy.loadtxt(zfile_name)
    valid = data[:,0]
    z_off = data[:,1]

    #
    # Load the locations of the beads.
    #
    #   This is a text file the contains the locations of the beads that 
    #   will be used to construct the PSF. Each line is a space separated 
    #   x, y pair of bead locations (in pixels).
    #
    #   One way to create this file is to look at the bead movie with
    #   visualizer.py and record the center positions of several beads.
    #
    data = numpy.loadtxt(beads_file, ndmin = 2)
    bead_x = data[:,0]
    bead_y = data[:,1]

    #
    # Go through the frames and the bead images to the average psf. Z 
    # positions are rounded to the nearest 50nm. You might need to 
    # adjust z_range depending on your experiment.
    #
    z_mid = int(z_range/z_step)
    max_z = 2 * z_mid + 1
    average_psf = numpy.zeros((max_z,4*aoi_size,4*aoi_size))
    totals = numpy.zeros(max_z)
    [dax_x, dax_y, dax_l] = movie_data.filmSize()
    for curf in range(dax_l):

        if ((curf%50)==0):
            print("Processing frame:", curf)

        if (abs(valid[curf]) < 1.0e-6):
            #    print "skipping", valid[curf]
            continue

        # Use bead localization to calculate spline.
        image = movie_data.loadAFrame(curf).astype(numpy.float64)

        # Get frame z and check that it is in range.
        zf = z_off[curf]
        zi = int(round(zf/z_step)) + z_mid
        if (zi > -1) and (zi < max_z):

            for i in range(bead_x.size):

                xf = bead_x[i]
                yf = bead_y[i]
                xi = int(xf)
                yi = int(yf)

                # Get localization image.
                mat = image[xi-aoi_size:xi+aoi_size,
                            yi-aoi_size:yi+aoi_size]
                
                # Zoom in by 2x.
                psf = scipy.ndimage.interpolation.zoom(mat, 2.0)

                # Re-center image.
                psf = scipy.ndimage.interpolation.shift(psf, (-2.0*(xf-xi), -2.0*(yf-yi)), mode='nearest')

                # Add to average psf accumulator.
                average_psf[zi,:,:] += psf
                totals[zi] += 1

    # Force PSF to be zero (on average) at the boundaries.
    for i in range(max_z):
        edge = numpy.concatenate((average_psf[i,0,:],
                                  average_psf[i,-1,:],
                                  average_psf[i,:,0],
                                  average_psf[i,:,-1]))
        average_psf[i,:,:] -= numpy.mean(edge)

    # Normalize PSF.
    for i in range(max_z):
        if (totals[i] > 0.0):
            average_psf[i,:,:] = average_psf[i,:,:]/numpy.sum(numpy.abs(average_psf[i,:,:]))

    average_psf = average_psf/numpy.max(average_psf)
    
    # Save PSF (in image form).
    if True:
        import storm_analysis.sa_library.daxwriter as daxwriter
        dxw = daxwriter.DaxWriter("psf_beads.dax", average_psf.shape[1], average_psf.shape[2])
        for i in range(max_z):
            #print i, numpy.max(average_psf[i,:,:])
            dxw.addFrame(1000.0 * average_psf[i,:,:] + 100)
        dxw.close()

    # Save PSF. 
    cur_z = -z_range
    z_vals = []
    for i in range(max_z):
        z_vals.append(cur_z)
        cur_z += z_step

    dict = {"psf" : average_psf,
            "pixel_size" : 0.080, # 1/2 the camera pixel size in nm.
            "type" : "3D",
            "zmin" : -z_range,
            "zmax" : z_range,
            "zvals" : z_vals}

    pickle.dump(dict, open(psf_name, "w"))


if (__name__ == "__main__"):

    import argparse
    
    parser = argparse.ArgumentParser(description = 'Measure PSF given a movie, a beads.txt file and a z_offset file')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--zoffset', dest='zoffset', type=str, required=True,
                        help = "A text file with two space separated numbers on each line, the first is 1 of the frame is valid, 0 otherwise and the second is the z offset of the frame relative to the focal plane in nanometers.")
    parser.add_argument('--beads', dest='beads', type=str, required=True,
                        help = "A text file with two space separated numbers on each line, the first is a bead X position in pixels and the second is a bead Y position")
    parser.add_argument('--psf', dest='psf', type=str, required=True,
                        help = "The name of the numpy format file to save the estimated PSF in.")
    parser.add_argument('--aoi_size', dest='aoi_size', type=int, required=False, default=12,
                        help = "The size of the area of interest around the bead in pixels. The default is 12.")
    parser.add_argument('--zrange', dest='zrange', type=float, required=False, default=750,
                        help = "The z range in nanometers. The PSF will be estimated from -zrange to +zrange. The default is 750nm.")
    parser.add_argument('--zstep', dest='zstep', type=float, required=False, default=50,
                        help = "The z step size in nanometers. The default is 50nm.")

    args = parser.parse_args()

    measurePSFBeads(args.movie, args.zoffset, args.beads, args.psf, aoi_size = args.aoi_size, z_range = args.zrange, z_step = args.zstep)
    
