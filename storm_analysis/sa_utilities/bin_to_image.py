#!/usr/bin/env python
"""
Functions to create images from a Insight3 localization binary file.

Note: This is deprecated as you should have HDF5 localization files.

Hazen 07/17
"""
import numpy
import sys
import tifffile

import storm_analysis.sa_library.grid_c as gridC
import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.i3togrid as i3togrid
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.simulator.draw_gaussians_c as dg


def render2DImage(i3_reader, shape, category = None, offsets = None, scale = 2, sigma = None):
    """
    Create a grayscale image from a Insight3 format binary data.

    i3_reader - A readinsight3.I3Reader object.
    shape - The shape of the output image. This will be multiplied by scale.
    category - Filter for localizations of this category. The default is all categories.
    offsets - X,Y offset of the image origin.  The default is 0.0.
    scale - The 'zoom' level of the output image, i.e. if the original STORM movie was
            256x256 and scale = 2 then the output image will be 512x512.
    sigma - The sigma to use when rendering gaussians (pixels). If this is None then
            the image will be a histogram.
    """
    image = numpy.zeros((shape[0]*scale, shape[1]*scale))

    # Make sure we are starting at the beginning.
    i3_reader.resetFp()

    # Load the first block of data.
    i3_data = i3_reader.nextBlock()
    while (i3_data is not False):
        sys.stdout.write(".")
        sys.stdout.flush()

        # Filter by category, if requested.
        if category is not None:
            i3_data = i3dtype.maskData(i3_data, (i3_data['c'] == category))

        # Adjust by offsets, if specified. Note, not adjusted for scale.
        if offsets is not None:
            i3_data['xc'] -= offsets[0]
            i3_data['yc'] -= offsets[1]

        # Adjust x,y by scale.
        xc = i3_data['xc']*scale
        yc = i3_data['yc']*scale

        # Histogram.
        if sigma is None:
            gridC.grid2D(numpy.round(xc),
                         numpy.round(yc),
                         image)
        # Gaussians.
        else:
            dg.drawGaussiansXYOnImage(image, xc, yc, sigma = sigma)

        # Load next block of data.
        i3_data = i3_reader.nextBlock()

    print()
    return image


def render2DImageFromFile(i3_filename, shape = None, category = None, offsets = None, scale = 2, sigma = None):
    """
    Wraps render2DIMage() to make it easier to create an 
    image from a Insight3 format binary file.
    """
    i3_reader = readinsight3.I3Reader(i3_filename)

    # If not specified, figure out shape.
    if shape is None:
        
        # Load the first block of data (which might be the whole file).
        i3_reader.resetFp()
        i3_data = i3_reader.nextBlock()

        # Try and figure out the original movie size.
        [x_size, y_size, temp] = i3togrid.getFilmSize(i3_filename, i3_data)

        shape = (x_size, y_size)

    return render2DImage(i3_reader,
                         shape,
                         category = category,
                         offsets = offsets,
                         scale = scale,
                         sigma = sigma)


def render3DImage(i3_reader, shape, category = None, offsets = None, scale = 2, sigma = None, z_edges = None):
    """
    Create a stack of grayscale images from a Insight3 format binary data.

    i3_reader - A readinsight3.I3Reader object.
    shape - The shape of the output image. This will be multiplied by scale.
    category - Filter for localizations of this category. The default is all categories.
    offsets - X,Y offset of the image origin.  The default is 0.0.
    scale - The 'zoom' level of the output image, i.e. if the original STORM movie was
            256x256 and scale = 2 then the output image will be 512x512.
    sigma - The sigma to use when rendering gaussians (pixels). If this is None then
            the image will be a histogram.
    z_edges - A list of z values specifying the z range for each image. This should be
            in nanometers.
    """
    num_z = len(z_edges)-1
    images = []
    for i in range(num_z):
        images.append(numpy.zeros((shape[0]*scale, shape[1]*scale)))

    # Make sure we are starting at the beginning.
    i3_reader.resetFp()
    
    # Load the first block of data.
    i3_data = i3_reader.nextBlock()
    while (i3_data is not False):
        sys.stdout.write(".")
        sys.stdout.flush()

        # Filter by category, if requested.
        if category is not None:
            i3_data = i3dtype.maskData(i3_data, (i3_data['c'] == category))

        # Adjust by offsets, if specified. Note, not adjusted for scale.
        if offsets is not None:
            i3_data['xc'] -= offsets[0]
            i3_data['yc'] -= offsets[1]

        # Adjust x,y by scale.
        xc = i3_data['xc']*scale
        yc = i3_data['yc']*scale

        # Iterate through z ranges.
        for i in range(num_z):
            z_mask = (i3_data['zc'] > z_edges[i]) & (i3_data['zc'] < z_edges[i+1])

            xc_z = xc[z_mask]
            yc_z = yc[z_mask]

            # Histogram.
            if sigma is None:
                gridC.grid2D(numpy.round(xc_z),
                             numpy.round(yc_z),
                             images[i])
            # Gaussians.
            else:
                dg.drawGaussiansXYOnImage(images[i], xc_z, yc_z, sigma = sigma)

        # Load next block of data.
        i3_data = i3_reader.nextBlock()

    print()
    return images

if (__name__ == "__main__"):

    import argparse
    
    parser = argparse.ArgumentParser(description = 'Create an image from an Insight3 format binary file.')

    parser.add_argument('--image', dest='image', type=str, required=True,
                        help = "The name of the output image (tiff format, 32 bit float).")
    parser.add_argument('--bin', dest='alist', type=str, required=True,
                        help = "The name of the localizations file. This is a binary file in Insight3 format.")
    parser.add_argument('--scale', dest='scale', type=int, required=False, default = 2,
                        help = "The 'zoom' of the output image (an integer).")
    parser.add_argument('--sigma', dest='sigma', type=float, required=False, default = 1.5,
                        help = "The sigma for gaussian render. Use 0.0 for a histogram.")    

    args = parser.parse_args()

    sigma = args.sigma
    if (sigma <= 0.0):
        sigma = None
        
    image = render2DImageFromFile(args.alist, scale = args.scale, sigma = sigma)

    tifffile.imsave(args.image, image.astype(numpy.float32))
