#!/usr/bin/env python
"""
Make images of the clusters.

Hazen 08/18
"""
import numpy
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import randomcolor

import storm_analysis.simulator.draw_gaussians_c as dg

import storm_analysis.dbscan.clusters_sa_h5py as clSAH5Py


def clusterImages(h5_name, min_size, image_max, scale = 4, show_unclustered = True):
    """
    Creates a RGB and an intensity array of the clusters with more than
    min_size elements.

    h5_name - The name of the HDF5 file with clustering information.
    min_size - The minimum size cluster to color.
    image_max - Maximum value for color normalization.
    scale - (int) Image rendering scale, default is 4.
    show_unclustered - Include unclustered localizations (in white).
    """
    rand_color = randomcolor.RandomColor()
    
    with clSAH5Py.SAH5Clusters(h5_name) as cl_h5:
        [movie_x, movie_y] = cl_h5.getMovieInformation()[:2]

        image_size = [movie_y * scale, movie_x * scale]

        red_image = numpy.zeros(image_size)
        grn_image = numpy.zeros(image_size)
        blu_image = numpy.zeros(image_size)
        sum_image = numpy.zeros(image_size)
    
        num_clusters = 0
        for index, cluster in cl_h5.clustersIterator(skip_unclustered = False, fields = ["x", "y"]):
            if ((index%200)==0):
                print("Processing cluster", index)

            x = scale * cluster['x']
            y = scale * cluster['y']

            if (index == 0) or (x.size < min_size):

                # Skip to next cluster if this cluster is too small and we are
                # not supposed to draw the unclustered localizations.
                if not show_unclustered:
                    continue

                else:
                    r = 1.0
                    g = 1.0
                    b = 1.0
            else:
                num_clusters += 1
                color = rand_color.generate()
                r = int(color[0][1:3], 16)/255.0
                g = int(color[0][3:5], 16)/255.0        
                b = int(color[0][5:7], 16)/255.0

            # FIXME: This would draw the image clusters as Gaussians. Not sure
            #        how well it actually works. Sigma should also probably be
            #        a user argument.
            #
            if False:
                temp = numpy.zeros(image_size)
                dg.drawGaussiansXYOnImage(temp, y, x, sigma = 1.5)

                red_image += r * temp
                grn_image += g * temp
                blu_image += b * temp
                sum_image += temp
        
            else:
                for i in range(x.size):
                    yi = int(x[i])
                    xi = int(y[i])
                    if (xi >= 0) and (xi < image_size[0]) and (yi >= 0) and (yi < image_size[1]):
                        red_image[xi,yi] += r
                        grn_image[xi,yi] += g
                        blu_image[xi,yi] += b
                        sum_image[xi,yi] += 1
        
        # Some wacky normalization scheme..
        mask = (sum_image > image_max)

        red_image[mask] = (red_image[mask]/sum_image[mask]) * image_max
        grn_image[mask] = (grn_image[mask]/sum_image[mask]) * image_max
        blu_image[mask] = (blu_image[mask]/sum_image[mask]) * image_max
        sum_image[mask] = image_max

        rgb_image = numpy.zeros((image_size[0], image_size[1], 3))
        rgb_image[:,:,0] = red_image
        rgb_image[:,:,1] = grn_image
        rgb_image[:,:,2] = blu_image

        rgb_image = (255.0/image_max)*rgb_image
        sum_image = (255.0/image_max)*sum_image

        rgb_image = rgb_image.astype(numpy.uint8)
        sum_image = sum_image.astype(numpy.uint8)

    return [rgb_image, sum_image, num_clusters]


def clusterImagesPNG(h5_name, png_base, title, min_size, image_max, scale):
    """
    Saves the results from clusterImages as a PNG file.

    FIXME: This was a convenience function that I used a lot in one paper
           but now it seems pretty specialized. Not sure it is worth 
           keeping in this project.
    """
    [rgb_image, sum_image, num_clusters] = clusterImages(h5_name, min_size, image_max, scale = scale)

    # Save RGB image.
    img1 = Image.fromarray(rgb_image, "RGB")
    try:
        draw1 = ImageDraw.Draw(img1)
        font1 = ImageFont.truetype("FreeMono.ttf", 24)
        draw1.text((2,2), title, (255,255,255), font = font1)
    except IOError:
        print("Text drawing disabled, true type font file may be missing?")
    
    img1.save(png_base + "_01.png")

    # Save grayscale image.
    img2 = Image.fromarray(sum_image, "L")
    try:
        img2 = img2.convert("RGB")
        draw2 = ImageDraw.Draw(img2)
        font2 = ImageFont.truetype("FreeMono.ttf", 24)
        draw2.text((2,2), title, (255,255,255), font = font2)
    except IOError:
        print("Text drawing disabled, true type font file may be missing?")
    
    img2.save(png_base + "_02.png")


    
if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Make images of the identified clusters in a cell')

    parser.add_argument('--bin', dest='clsah5', type=str, required=True,
                        help = "The name of the clustered HDF5 format localization file.")
    parser.add_argument('--image', dest='image', type=str, required=True,
                        help = "Root name of image file. Two images are generated, one with the clusters colored and one without.")
    parser.add_argument('--title', dest='title', type=str, required=True,
                        help = "Title to add to the output images.")
    parser.add_argument('--min_size', dest='min_size', type=int, required=True,
                        help = "The minimum cluster size to render.")
    parser.add_argument('--im_max', dest='image_max', type=int, required=True,
                        help = "Image maximum, a good value for this is (usually) something in the range 6-20.")
    parser.add_argument('--scale', dest='scale', type=int, required=False, default = 4,
                        help = "The 'zoom' of the output image (an integer).")
    
    args = parser.parse_args()
    
    clusterImages(args.clsah5, args.image, args.title, args.min_size, args.image_max, args.scale)

