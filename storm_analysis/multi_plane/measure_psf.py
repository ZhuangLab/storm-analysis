#!/usr/bin/env python
"""
Measure the PSF given the averaged z-stacks created using  
multi_plane.psf_zstack and a text file containing the z offsets 
of the images in the z-stack as created by spliner.offset_to_z.

This is pretty similar to spliner, except that we do not normalize
each Z section such that its (absolute) value sums to 1.0. Instead
each Z section is normalized by the number of events / observations.

Note: For multi-color imaging you probably want to use normalize
      True to normalize each PSF independently.

FIXME: Need to be able to adjust for different starting frames, 
       as the cameras will not necessarily all start at the 
       same time.

Hazen 05/17
"""

import pickle
import numpy
import os
import tifffile

import storm_analysis.spliner.measure_psf_utils as measurePSFUtils

def measurePSF(zstack_name, zfile_name, psf_name, pixel_size = 0.1, refine = False, z_range = 0.75, z_step = 0.050, normalize = False):
    """
    zstack_name - The name of the file containing the z-stacks.
    zfile_name - The text file containing the z offsets (in microns) for each frame.
    psf_name - The name of the file to save the measured PSF in (as a pickled Python dictionary).
    pixel_size - The pixel size in microns.
    refine - Align the measured PSF for each bead to the average PSF.
    z_range - The range the PSF should cover in microns.
    z_step - The z step size of the PSF.
    normalize - If true, normalize the PSF to unit height.
    """

    # Create z scaling object.
    z_sclr = measurePSFUtils.ZScaler(z_range, z_step)
    max_z = z_sclr.getMaxZ()

    # Load z-stacks.
    zstacks = numpy.load(zstack_name)
    x_size = zstacks[0].shape[0]
    y_size = zstacks[0].shape[1]

    # Load z-offsets.
    z_offset_data = numpy.loadtxt(zfile_name, ndmin = 2)
    is_valid = z_offset_data[:,0]
    z_offsets = z_offset_data[:,1]

    # Check if the z-stack has fewer frames than the offset file.
    n_frames = z_offsets.size
    if (n_frames > zstacks[0].shape[2]):
        print("Warning! z stack has", n_frames - zstacks[0].shape[2], "fewer frames than the z offset file.")
        n_frames = zstacks[0].shape[2]
    
    # Average stacks in z.
    psfs = []
    for i in range(len(zstacks)):
        psf = numpy.zeros((max_z, x_size, y_size))
        totals = numpy.zeros(max_z, dtype = numpy.int)

        for j in range(n_frames):

            # 0.0 = invalid, 1.0 = valid.
            if (is_valid[j] > 1.0e-3):
                zi = z_sclr.convert(z_offsets[j])
                if z_sclr.inRange(zi):
                    psf[zi,:,:] += zstacks[i][:,:,j]
                    totals[zi] += 1

        # Check that we got at least one valid measurement. Totals
        # is expected to be the same for each PSF.
        if (i == 0):
            assert (numpy.max(totals) > 0)

        # Normalize each PSF z plane by the total counts in the plane.
        for j in range(max_z):
            if (i == 0):
                print("z plane {0:0d} has {1:0d} samples".format(j, totals[j]))
            if (totals[j] > 0):
                psf[j,:,:] = psf[j,:,:]/float(totals[j])

        # Append to list of PSFs.
        psfs.append(psf)

    # Align the PSFs to each other. This should hopefully correct for
    # any small errors in the input locations, and also for fields of
    # view that are not completely flat.
    #
    if refine:
        print("Refining PSF alignment.")
        [average_psf, i_score] = measurePSFUtils.alignPSFs(psfs)
    else:
        average_psf = measurePSFUtils.averagePSF(psfs)        

    # Normalize to unity height, if requested.
    if normalize and (numpy.amax(average_psf) > 0.0):
        print("Normalizing PSF.")
        average_psf = average_psf/numpy.amax(average_psf)

    if not (numpy.amax(average_psf) > 0.0):
        print("Warning! Measured PSF maxima is zero or negative!")

    # Save PSF.
    #
    #  At least for now the PSFs use nanometers, not microns.
    #
    z_range = z_range * 1.0e+3
    z_step = z_step * 1.0e+3        

    cur_z = -z_range
    z_vals = []
    for i in range(max_z):
        z_vals.append(cur_z)
        cur_z += z_step

    psf_dict = {"maximum" : numpy.amax(average_psf),
                "psf" : average_psf,
                "pixel_size" : pixel_size,
                "type" : "3D",
                "version" : 2.0,
                "zmin" : -z_range,
                "zmax" : z_range,
                "zvals" : z_vals}

    with open(psf_name, 'wb') as fp:
        pickle.dump(psf_dict, fp)

    # Save (normalized) z_stack as tif for inspection purposes.
    if (numpy.amax(average_psf) > 0.0):
        average_psf = average_psf/numpy.amax(average_psf)
    average_psf = average_psf.astype(numpy.float32)
    with tifffile.TiffWriter(os.path.splitext(psf_name)[0] + ".tif") as tf:
        for i in range(max_z):
            tf.save(average_psf[i,:,:])


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Measure PSF from an average z-stack and a z offset file.')

    parser.add_argument('--zstack', dest='zstack', type=str, required=True,
                        help = "The name of the numpy file containing the averaged z-stack.")
    parser.add_argument('--zoffsets', dest='zoffsets', type=str, required=True,
                        help = "The name of the text file containing the per-frame z offsets (in microns).")
    parser.add_argument('--psf_name', dest='psf_name', type=str, required=True,
                        help = "The name of the file for saving the measured PSF.")
    parser.add_argument('--pixel_size', dest='pixel_size', type=float, required=False, default=100.0,
                        help = "The pixel size in nanometers. The default is 100nm.")
    parser.add_argument('--refine', dest='refine', action='store_true', default=False,
                        help = "Refine PSF locations relative to each other, default is False")
    parser.add_argument('--z_range', dest='z_range', type=float, required=False, default=0.75,
                        help = "The z range (+-) in microns, default is +-0.75um.")
    parser.add_argument('--z_step', dest='z_step', type=float, required=False, default=0.05,
                        help = "The z step size in microns. The default is 0.05um.")
    parser.add_argument('--normalize', dest='norm', action='store_true', default=False,
                        help = "Normalize PSFs to unit height, default is False")

    args = parser.parse_args()
    
    measurePSF(args.zstack,
               args.zoffsets,
               args.psf_name,
               pixel_size = args.pixel_size * 1.0e-3,
               refine = args.refine,
               z_range = args.z_range,
               z_step = args.z_step,
               normalize = args.norm)
