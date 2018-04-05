#!/usr/bin/env python
"""
Make pictures of a PSF.

Hazen 04/18
"""
import pickle
import matplotlib
import matplotlib.pyplot as pyplot


def psfImages(psf_filename, verbose = True, sx = 12, sy = 4):
    
    with open(psf_filename, 'rb') as fp:
        psf_data = pickle.load(fp)       

    psf = psf_data["psf"]

    if verbose:
        print("PSF shape:", psf.shape)
        print("pixel size: {0:.3f}um".format(psf_data["pixel_size"]))
        print("zmin, zmax: {0:.1f}nm, {1:.1f}nm".format(psf_data["zmin"], psf_data["zmax"]))

    mid_xy = int(psf.shape[1]/2)
    mid_z = int(psf.shape[0]/2)
    xy_max = psf_data["pixel_size"] * psf.shape[1]
    xy_min = 0.0
    z_min = psf_data["zmin"] * 1.0e-3
    z_max = psf_data["zmax"] * 1.0e-3

    fig = pyplot.figure(figsize = (12,4))
    ax1 = fig.add_subplot(1,3,1)
    ax1.imshow(psf[mid_z,:,:], 
               interpolation = 'none', 
               extent = [xy_min, xy_max, xy_min, xy_max],
               cmap = "gray")
    ax1.set_title("PSF XY slice")

    ax2 = fig.add_subplot(1,3,2)
    ax2.imshow(psf[:,mid_xy,:],
               interpolation = 'none', 
               extent = [xy_min, xy_max, z_min, z_max],
               cmap = "gray")
    ax2.set_title("PSF YZ slice")

    ax3 = fig.add_subplot(1,3,3)
    ax3.imshow(psf[:,:,mid_xy], 
               interpolation = 'none', 
               extent = [xy_min, xy_max, z_min, z_max],
               cmap = "gray")
    ax3.set_title("PSF XZ slice")

    pyplot.show()

    if verbose:
        print("Plots are in microns")


if (__name__ == "__main__"):
    
    import argparse

    parser = argparse.ArgumentParser(description = 'Make images of a PSF.')

    parser.add_argument('--psf', dest='psf', type=str, required=True,
                        help = "The name of the PSF format file.")

    args = parser.parse_args()

    psfImages(args.psf)
