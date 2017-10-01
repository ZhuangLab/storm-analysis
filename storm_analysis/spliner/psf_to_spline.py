#!/usr/bin/env python
"""
Takes a 3D numpy array (as created by measure_psf.py) and
outputs a 3D numpy array that can be used as a spline for
3D spline fitting.

Size is the size of the spline in pixels. This should be
a multiple of 2.

Hazen 01/14
"""

import os
import pickle
import numpy
import sys

import storm_analysis.spliner.spline1D as spline1D
import storm_analysis.spliner.spline2D as spline2D
import storm_analysis.spliner.spline3D as spline3D


def psfToSpline(psf_name, spline_name, s_size):

    # Load PSF.
    with open(psf_name, 'rb') as fp:
        psf_data = pickle.load(fp)        
    np_psf = psf_data["psf"]
    
    spline = False
    start = np_psf.shape[1]/2.0 - s_size - 0.5

    # 2D spline
    if (len(np_psf.shape) == 2):
        print("Generating 2D spline.")
        s_size = 2*s_size

        np_spline = numpy.zeros((s_size, s_size))
        #np_psf = np_psf/numpy.max(np_psf)
        xy_spline = spline2D.Spline2D(np_psf)
        
        x = start
        for i in range(s_size):
            y = start
            for j in range(s_size):
                np_spline[j,i] = xy_spline.f(y,x)
            
                y += 1.0
            x += 1.0

        print("Calculating spline coefficients.")
        spline = spline2D.Spline2D(np_spline)

        if True:
            import tifffile
            tiff_name = os.path.splitext(spline_name)[0] + "_sp.tif"
            tifffile.imsave(tiff_name, np_spline.astype(numpy.float32))


    # 3D spline
    else:
        print("Generating 3D spline.")
        s_size = 2*s_size

        np_spline = numpy.zeros((s_size, s_size, s_size))
        xy_splines = []

        print("Generating XY splines.")
        for i in range(np_psf.shape[0]):
            xy_splines.append(spline2D.Spline2D(np_psf[i,:,:]))

        print("Generating fitting spline.")
        x = start
        for i in range(s_size):
            y = start
            for j in range(s_size):

                zvals = numpy.zeros(np_psf.shape[0])
                for k in range(np_psf.shape[0]):
                    zvals[k] = xy_splines[k].f(y,x)
                    z_spline = spline1D.Spline1D(zvals)

                max_z = float(np_psf.shape[0]) - 1.0
                inc = max_z/(float(s_size)-1.0)
                for k in range(s_size):
                    z = float(k)*inc
                    if (z > max_z):
                        z = max_z
                    np_spline[k,j,i] = z_spline.f(z)

                y += 1.0
            x += 1.0

        print("Calculating spline coefficients.")
        spline = spline3D.Spline3D(np_spline, verbose = True)

        if True:
            import tifffile
            tiff_name = os.path.splitext(spline_name)[0] + "_sp.tif"
            with tifffile.TiffWriter(tiff_name) as tf:
                for i in range(s_size):
                    tf.save(np_spline[i,:,:].astype(numpy.float32))

    del psf_data["psf"]
    psf_data["spline"] = np_spline
    psf_data["coeff"] = spline.getCoeff()
    psf_data["psf_name"] = psf_name
    with open(spline_name, 'wb') as fp:
        pickle.dump(psf_data, fp)


if (__name__ == "__main__"):

    import argparse
    
    parser = argparse.ArgumentParser(description = 'Convert PSF into a spline file.')

    parser.add_argument('--psf', dest='psf', type=str, required=True,
                        help = "The name of the numpy format file containing the estimated PSF.")
    parser.add_argument('--spline', dest='spline', type=str, required=True,
                        help = "The name of the numpy format file to save the spline in.")
    parser.add_argument('--spline_size', dest='spline_size', type=int, required=True,
                        help = "The size of the spline in pixels in X, Y")
    
    args = parser.parse_args()
    
    psfToSpline(args.psf, args.spline, args.spline_size)


