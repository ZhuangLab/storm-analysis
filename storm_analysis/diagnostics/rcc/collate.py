#!/usr/bin/env python
"""
Check drift correction results.

Hazen 01/18
"""
import glob
import numpy
import tifffile

import storm_analysis.sa_utilities.hdf5_to_image as hdf5ToImage

import storm_analysis.diagnostics.rcc.settings as settings


def collate():
    # We assume that there are two directories, the first with xy correction
    # only and the second with xyz correction.
    #

    # Make 2D images of the drift corrected data.
    if False:
        for adir in ["test_01/", "test_02/"]:
            im = hdf5ToImage.render2DImage(adir + "test.hdf5", sigma = 0.5)
            tifffile.imsave(adir + "test_im_2d.tif", im.astype(numpy.float32))

    # Make 3D images of the drift corrected data.
    if False:
        z_edges = numpy.arange(-0.5, 0.55, 0.1)
        for adir in ["test_01/", "test_02/"]:
            images = hdf5ToImage.render3DImage(adir + "test.hdf5", z_edges, sigma = 0.5)
            with tifffile.TiffWriter(adir + "test_im_3d.tif") as tf:
                for elt in images:
                    tf.save(elt.astype(numpy.float32))

    # Measure error in drift measurements.
    if True:
        print("Drift correction RMS error (nanometers):")
        for i, elt in enumerate(["drift_xy.txt", "drift_xyz.txt"]):

            print(" ", elt)
    
            ref = numpy.loadtxt(elt)
            exp = numpy.loadtxt("test_{0:02}/test_drift.txt".format(i+1))
        
            max_len = exp.shape[0]
            for j, elt in enumerate(["X", "Y", "Z"]):
                refv = ref[:max_len,j]
                expv = exp[:,j+1]

                # Correct for DC offset.
                expv += numpy.mean(refv - expv)

                # Calculate error (in nanometers).
                if (elt == "Z"):
                    print("  ", elt, "{0:.3f}".format(numpy.std(refv - expv) * 1.0e+3))
                else:
                    print("  ", elt, "{0:.3f}".format(numpy.std(refv - expv) * settings.pixel_size))
            print()


if (__name__ == "__main__"):
    collate()
    
        
        
