#!/usr/bin/env python
"""
Measure the PSF given the averaged, 2x up-sampled z-stacks
created using multi_plane.psf_zstack and a text file containing
the z offsets of the images in the z-stack as created by
spliner.offset_to_z.py.

The output is pretty similar to spliner.measure_psf, except that 
this does not normalize to unity by default as differences in 
intensity between the different image planes are expected.

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


def measurePSF(zstack_name, zfile_name, psf_name, z_range = 750.0, z_step = 50.0, normalize = False):

    # Load z-stack.
    zstack = numpy.load(zstack_name)
    x_size = zstack.shape[0]
    y_size = zstack.shape[1]

    # Load z-offsets.
    z_offset_data = numpy.loadtxt(zfile_name, ndmin = 2)
    is_valid = z_offset_data[:,0]
    z_offsets = z_offset_data[:,1]

    # Check if the z-stack has fewer frames than the offset file.
    n_frames = z_offsets.size
    if (n_frames > zstack.shape[2]):
        print("Warning! z stack has", n_frames - zstack.shape[2], "fewer frames than the z offset file.")
        n_frames = zstack.shape[2]
    
    # Average stack in z.
    z_mid = int(z_range/z_step)
    max_z = 2 * z_mid + 1

    average_psf = numpy.zeros((max_z, x_size, y_size))
    totals = numpy.zeros(max_z)
    for i in range(n_frames):

        # 0.0 = in valid, 1.0 = valid.
        if (is_valid[i] > 1.0e-3):
            zi = int(round(z_offsets[i]/z_step) + z_mid)
            if (zi > -1) and (zi < max_z):
                average_psf[zi,:,:] += zstack[:,:,i]
                totals[zi] += 1

    # Normalize each z plane by the total counts in the plane.
    for i in range(max_z):
        print(i, totals[i])
        if (totals[i] > 0):
            average_psf[i,:,:] = average_psf[i,:,:]/totals[i]

    # Normalize to unity height, if requested.
    if normalize:
        print("Normalizing PSF.")
        average_psf = average_psf/numpy.amax(average_psf)
    
    # Save PSF.
    cur_z = -z_range
    z_vals = []
    for i in range(max_z):
        z_vals.append(cur_z)
        cur_z += z_step

    psf_dict = {"maximum" : numpy.amax(average_psf),
                "psf" : average_psf,
                "type" : "3D",
                "zmin" : -z_range,
                "zmax" : z_range,
                "zvals" : z_vals}

    with open(psf_name, 'wb') as fp:
        pickle.dump(psf_dict, fp)

    # Save (normalized) z_stack as tif for inspection purposes.
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
                        help = "The name of the text file containing the per-frame z offsets (in nm).")
    parser.add_argument('--psf_name', dest='psf_name', type=str, required=True,
                        help = "The name of the file for saving the measured PSF.")
    parser.add_argument('--z_range', dest='z_range', type=float, required=False, default=750.0,
                        help = "The z range (+-) in nm, default is +-750nm.")
    parser.add_argument('--z_step', dest='z_step', type=float, required=False, default=50.0,
                        help = "The z step size in nm. The default is 50nm.")
    parser.add_argument('--normalize', dest='norm', type=bool, required=False, default=False,
                        help = "Normalize PSF height to unity.")

    args = parser.parse_args()
    
    measurePSF(args.zstack,
               args.zoffsets,
               args.psf_name,
               z_range = args.z_range,
               z_step = args.z_step,
               normalize = args.norm)
