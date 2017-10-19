#!/usr/bin/env python
"""
Converts a PSF measured by Spliner into one that is
compatible with PSF FFT.

Hazen 10/17
"""
import pickle
import math
import numpy
import tifffile

import storm_analysis.sa_library.rebin as rebin

def downsamplePSF(spliner_psf_filename, psf_filename, pixel_size):
    """
    pixel_size - The final pixel_size (in microns). This does not actually 
                 change anything, it just makes it easy to change the pixel 
                 size to match that of the analysis.xml file so that PSF 
                 FFT does not complain.
    """
    with open(spliner_psf_filename, 'rb') as fp:
        spliner_psf_data = pickle.load(fp)

    spliner_psf_data["pixel_size"] = pixel_size

    psf = spliner_psf_data["psf"]

    psf = rebin.downSample(psf,
                           psf.shape[0],
                           int(psf.shape[1]/2),
                           int(psf.shape[2]/2))

    spliner_psf_data["psf"] = psf

    with open(psf_filename, 'wb') as fp:
        pickle.dump(spliner_psf_data, fp)

    # Also save a .tif version.
    with tifffile.TiffWriter("psf.tif") as tf:
        for i in range(psf.shape[0]):
            tf.save(psf[i,:,:].astype(numpy.float32))   


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Downsample a Spliner PSF.')

    parser.add_argument('--spliner_psf', dest='spliner_psf', type=str, required=True,
                        help = "The name of the file containing the PSF as measured by Spliner.")
    parser.add_argument('--psf', dest='psf', type=str, required=True,
                        help = "The name of the file for saving the downsampled PSF.")
    parser.add_argument('--pixel-size', dest='pixel_size', type=float, required=True,
                        help = "The pixel size in nanometers.")

    args = parser.parse_args()

    downsamplePSF(args.spliner_psf, args.psf, args.pixel_size * 1.0e-3)
