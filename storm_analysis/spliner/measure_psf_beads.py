#!/usr/bin/env python
"""
Measure the 3D PSF a movie given the locations of the beads
of interest in the movie and the z-offset of each frame of
the movie. It is assumed that the drift over the time
course of the movie is neglible.

Depending on your setup you may need to change:
  1. The z range (z_range).
  2. The pixel size (pixel_size).
  3. The AOI size (aoi_size). This is less important as you
     will get to specify the final size to use when you use
     psf_to_spline.py to create the spline to use for fitting.

Hazen 1/16
"""
import os
import pickle
import numpy
import scipy
import scipy.ndimage
import sys
import tifffile

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.parameters as params
import storm_analysis.spliner.measure_psf_utils as measurePSFUtils


def measurePSFBeads(movie_name, zfile_name, beads_file, psf_name, aoi_size = 12, pixel_size = 0.1, refine = False, z_range = 0.6, z_step = 0.05):
    """
    movie_name - The name of the movie, presumably a z stack for PSF measurement.
    zfile_name - The text file containing the z offsets (in microns) for each frame.
    beads_file - The text file containing the locations of the beads.
    psf_name - The name of the file to save the measured PSF in (as a pickled Python dictionary).
    aoi_size - The AOI of interest in pixels. The final AOI is 2x this number.
    pixel_size - The pixel size in microns.
    refine - Align the measured PSF for each bead to the average PSF.
    z_range - The range the PSF should cover in microns.
    z_step - The z step size of the PSF.
    """
    # Load the z-offset information for the dax file.
    #
    #   This is a text file with one line per frame that contains the 
    #   z-offset (in microns) for that frame. Each line is a space separated
    #   valid, z_pos pair. If valid if 0 the frame will be ignored,
    #   otherwise it will be used.
    #
    z_offsets = numpy.loadtxt(zfile_name)

    # Create array specifying what frame corresponds to what
    # Z slice in the PSF.
    #
    z_index = measurePSFUtils.makeZIndexArray(z_offsets, z_range, z_step)

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
    bead_x = data[:,1] + 1
    bead_y = data[:,0] + 1    

    # Create a reader of the movie.
    #
    #   We assume that the bead stack was measured with a camera
    #   that does not have a large pixel to pixel variation in
    #   gain and offset. The offset and magnitude are not that
    #   important at we will estimate and subtract the offset
    #   and normalize 1.
    #

    # Movie (frame) reader.
    frame_reader = analysisIO.FrameReaderStd(movie_file = movie_name,
                                             camera_gain = 1.0,
                                             camera_offset = 0.0)
    
    # Measure PSFs for each bead.
    #
    total_samples = None
    psfs = []
    for i in range(bead_x.size):
        [psf, samples] = measurePSFUtils.measureSinglePSFBeads(frame_reader,
                                                               z_index,
                                                               aoi_size,
                                                               bead_x[i],
                                                               bead_y[i],
                                                               zoom = 1)

        # Verify that we have at least one sample per section, because if
        # we don't this almost surely means something is wrong.
        if (i == 0):
            for j in range(samples.size):
                assert(samples[j] > 0), "No data for PSF z section " + str(j)
        
        # Normalize by the number of sample per z section.
        #for j in range(samples.size):
        #    psf[j,:,:] = psf[j,:,:]/samples[j]

        # Keep track of total number of samples.
        if total_samples is None:
            total_samples = samples
        else:
            total_samples += samples

        psfs.append(psf)

    # Set the PSF to have zero average on the X/Y boundaries. We are
    # matching the behavior of spliner.measure_psf here.
    #
    sum_psf = measurePSFUtils.sumPSF(psfs)
    for i in range(sum_psf.shape[0]):
        mean_edge = measurePSFUtils.meanEdge(sum_psf[i,:,:])
        for j in range(len(psfs)):
            psfs[j][i,:,:] -= mean_edge/float(len(psfs))

    # Align the PSFs to each other. This should hopefully correct for
    # any small errors in the input locations, and also for fields of
    # view that are not completely flat.
    #
    if refine:
        print("Refining PSF alignment.")

        # Normalize each PSF by the number of z sections.
        for psf in psfs:
            for i in range(samples.size):
                psf[i,:,:] = psf[i,:,:]/samples[i]
            
        [average_psf, i_score] = measurePSFUtils.alignPSFs(psfs)
    else:
        average_psf = measurePSFUtils.averagePSF(psfs)

    # Normalize PSF.
    #
    #   This normalizes the PSF so that sum of the absolute values
    #   of each section is 1.0. This only makes sense if the AOI is
    #   large enough to capture all the photons, which might not be
    #   true. Not clear how important this is as Spliner will fit
    #   for the height anyway.
    #
    for i in range(average_psf.shape[0]):
        print("z plane {0:0d} has {1:0d} samples".format(i, total_samples[i]))
        
        section_sum = numpy.sum(numpy.abs(average_psf[i,:,:]))
        
        # Do we need this test? We already check that we have at
        # least one sample per section.
        if (section_sum > 0.0):
            average_psf[i,:,:] = average_psf[i,:,:]/section_sum

    # Normalize to unity maximum height.
    if (numpy.max(average_psf) > 0.0):
        average_psf = average_psf/numpy.max(average_psf)
    else:
        print("Warning! Measured PSF maxima is zero or negative!")
    
    # Save PSF (in image form).
    if True:
        tif_name = os.path.splitext(psf_name)[0]
        with tifffile.TiffWriter(tif_name + "_beads.tif") as tf:
            for i in range(average_psf.shape[0]):
                tf.save(average_psf[i,:,:].astype(numpy.float32))

    # Save PSF.
    #
    #   For historical reasons all the PSF z values are in nanometers.
    #   At some point this should be fixed.
    #
    z_range = 1.0e+3 * z_range
    z_step = 1.0e+3 * z_step
    
    cur_z = -z_range
    z_vals = []
    for i in range(average_psf.shape[0]):
        z_vals.append(cur_z)
        cur_z += z_step

    psf_dict = {"psf" : average_psf,
                "pixel_size" : pixel_size,
                "type" : "3D",
                "version" : 2.0,
                "zmin" : -z_range,
                "zmax" : z_range,
                "zvals" : z_vals}

    pickle.dump(psf_dict, open(psf_name, 'wb'))


if (__name__ == "__main__"):

    import argparse
    
    parser = argparse.ArgumentParser(description = 'Measure PSF given a movie, a beads.txt file and a z_offset file')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--zoffset', dest='zoffset', type=str, required=True,
                        help = "A text file with two space separated numbers on each line, the first is 1 of the frame is valid, 0 otherwise and the second is the z offset of the frame relative to the focal plane in microns.")
    parser.add_argument('--beads', dest='beads', type=str, required=True,
                        help = "A text file with two space separated numbers on each line, the first is a bead X position in pixels and the second is a bead Y position")
    parser.add_argument('--psf', dest='psf', type=str, required=True,
                        help = "The name of the numpy format file to save the estimated PSF in.")
    parser.add_argument('--aoi_size', dest='aoi_size', type=int, required=False, default=12,
                        help = "The size of the area of interest around the bead in pixels. The default is 12.")
    parser.add_argument('--pixel_size', dest='pixel_size', type=float, required=False, default=100.0,
                        help = "The movie pixel size in nanometers. The default is 100nm.")
    parser.add_argument('--refine', dest='refine', action='store_true', default=False)
    parser.add_argument('--zrange', dest='zrange', type=float, required=False, default=0.75,
                        help = "The z range in microns. The PSF will be estimated from -zrange to +zrange. The default is 0.75um.")
    parser.add_argument('--zstep', dest='zstep', type=float, required=False, default=0.05,
                        help = "The z step size in microns. The default is 0.05um.")
    
    args = parser.parse_args()

    measurePSFBeads(args.movie,
                    args.zoffset,
                    args.beads,
                    args.psf,
                    aoi_size = args.aoi_size,
                    pixel_size = args.pixel_size * 1.0e-3,
                    refine = args.refine,
                    z_range = args.zrange,
                    z_step = args.zstep)
    
