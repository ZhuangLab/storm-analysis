#!/usr/bin/env python
"""
Functions to create images from a HDF5 localization binary file. The
orientations of the output images are choosen such that when saved
using tifffile they match that of visualizer/visualizer.py.

Hazen 01/18
"""
import numpy
import sys

import storm_analysis.sa_library.grid_c as gridC
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.simulator.draw_gaussians_c as dg


def filterOffsetScale(locs, category, offsets, scale):
    """
    Filter by category, offset and scale.
    """
    # Filter by category, if requested.
    if (category is not None):
        mask = (locs["category"] == category)

        if (numpy.count_nonzero(mask) == 0):
            return None
        
        locs["x"] = locs["x"][mask]
        locs["y"] = locs["y"][mask]
        if "z" in locs:
            locs["z"] = locs["z"][mask]

    # Adjust by offsets, if specified. Note, not adjusted for scale.        
    if offsets is not None:
        locs["x"] += offsets[0]
        locs["y"] += offsets[1]

    # Adjust x,y by scale.
    locs["x"] = locs["x"] * scale
    locs["y"] = locs["y"] * scale

    return locs

        
def renderImage(image, x, y, sigma):
    """
    A helper function that does the actual rendering.
    """
    # Histogram.
    if sigma is None:
        gridC.grid2D(numpy.round(y),
                     numpy.round(x),
                     image)
    # Gaussians.
    else:
        dg.drawGaussiansXYOnImage(image, y, x, sigma = sigma)
        
        
def render2DImage(h5_name, category = None, offsets = None, scale = 2, sigma = None):
    """
    Create a grayscale image from a HDF5 format localizations file. This will use
    the tracks if available, otherwise it will use the localizations.

    h5_name - The name of the HDF5 file.
    category - Filter for localizations of this category. The default is all categories.
    offsets - List containing [X,Y] offset of the image origin.  The default is no offset.
    scale - The 'zoom' level of the output image, i.e. if the original STORM movie was
            256x256 and scale = 2 then the output image will be 512x512.
    sigma - The sigma to use when rendering gaussians (pixels). If this is None then
            the image will be a histogram.
    """
    with saH5Py.SAH5Py(h5_name) as h5:
        [movie_x, movie_y, movie_l, hash_value] = h5.getMovieInformation()

        if sigma is None:
            image = numpy.zeros((movie_y * scale, movie_x * scale), dtype = numpy.int32)
        else:
            image = numpy.zeros((movie_y * scale, movie_x * scale), dtype = numpy.float64)

        fields = ["x", "y"]
        if category is not None:
            fields.append("category")

        if h5.hasTracks():
            for locs in h5.tracksIterator(fields = fields):
                sys.stdout.write(".")
                sys.stdout.flush()
                locs = filterOffsetScale(locs, category, offsets, scale)
                if locs is not None:
                    renderImage(image, locs["x"], locs["y"], sigma)
            
        else:
            print("Tracks not found, using localizations.")
            for fnum, locs in h5.localizationsIterator(fields = fields):
                if ((fnum%2000)==0):
                    sys.stdout.write(".")
                    sys.stdout.flush()
                locs = filterOffsetScale(locs, category, offsets, scale)
                if locs is not None:
                    renderImage(image, locs["x"], locs["y"], sigma)

    print()
    return image


def render3DImage(h5_name, z_edges, category = None, offsets = None, scale = 2, sigma = None):
    """
    Create a stack of grayscale image from a HDF5 format localizations file. This will use
    the tracks if available, otherwise it will use the localizations.

    h5_name - The name of the HDF5 file.
    category - Filter for localizations of this category. The default is all categories.
    offsets - X,Y offset of the image origin.  The default is 0.0.
    scale - The 'zoom' level of the output image, i.e. if the original STORM movie was
            256x256 and scale = 2 then the output image will be 512x512.
    sigma - The sigma to use when rendering gaussians (pixels). If this is None then
            the image will be a histogram.
    z_edges - A list of z values specifying the z range for each image. This should be
            in microns.
    """
    with saH5Py.SAH5Py(h5_name) as h5:
        [movie_x, movie_y, movie_l, hash_value] = h5.getMovieInformation()

        num_z = len(z_edges)-1
        images = []
        for i in range(num_z):
            if sigma is None:
                images.append(numpy.zeros((movie_y * scale, movie_x * scale), dtype = numpy.int32))
            else:
                images.append(numpy.zeros((movie_y * scale, movie_x * scale), dtype = numpy.float64))

        fields = ["x", "y", "z"]
        if category is not None:
            fields.append("category")
            
        if h5.hasTracks():
            for locs in h5.tracksIterator(fields = fields):
                sys.stdout.write(".")
                sys.stdout.flush()
                locs = filterOffsetScale(locs, category, offsets, scale)
                
                if locs is not None:
                    
                    # Iterate through z ranges.
                    for i in range(num_z):
                        z_mask = (locs["z"] >= z_edges[i]) & (locs["z"] < z_edges[i+1])
                        renderImage(images[i],
                                    locs["x"][z_mask],
                                    locs["y"][z_mask],
                                    sigma)
            
        else:
            print("Tracks not found, using localizations.")
            for fnum, locs in h5.localizationsIterator(fields = fields):
                if ((fnum%2000)==0):
                    sys.stdout.write(".")
                    sys.stdout.flush()
                    
                locs = filterOffsetScale(locs, category, offsets, scale)    
                if locs is not None:
                    
                    # Iterate through z ranges.
                    for i in range(num_z):
                        z_mask = (locs["z"] >= z_edges[i]) & (locs["z"] < z_edges[i+1])
                        renderImage(images[i],
                                    locs["x"][z_mask],
                                    locs["y"][z_mask],
                                    sigma)

    print()
    return images


if (__name__ == "__main__"):

    import argparse
    import tifffile
    
    parser = argparse.ArgumentParser(description = 'Create a 2D image from an HDF5 format localization file.')

    parser.add_argument('--image', dest='image', type=str, required=True,
                        help = "The name of the output image (tiff format, 32 bit float).")
    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the localizations HDF5 file.")
    parser.add_argument('--scale', dest='scale', type=int, required=False, default = 2,
                        help = "The 'zoom' of the output image (an integer).")
    parser.add_argument('--sigma', dest='sigma', type=float, required=False, default = 1.5,
                        help = "The sigma for gaussian render. Use 0.0 for a histogram.")    

    args = parser.parse_args()

    sigma = args.sigma
    if (sigma <= 0.0):
        sigma = None
        
    image = render2DImage(args.hdf5, scale = args.scale, sigma = sigma)

    # Get image pixel size.
    pixel_size = None
    try:
        with saH5Py.SAH5Py(args.hdf5) as h5:
            pixel_size = h5.getPixelSize()/args.scale
    except saH5Py.SAH5PyException:
        print("No pixel size information available.")

    with tifffile.TiffWriter(open(args.image, "wb")) as tf:
        if pixel_size is not None:
            resolution = [1.0e+7/pixel_size, 1.0e+7/pixel_size, 3]
            tf.save(image.astype(numpy.float32), resolution=resolution)
        else:
            tf.save(image.astype(numpy.float32))

