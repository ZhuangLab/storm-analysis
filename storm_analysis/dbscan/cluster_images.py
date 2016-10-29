#!/usr/bin/python
#
# Make images of the identified clusters in a cell.
#
# Hazen 09/16
#

import numpy
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import randomcolor
import sys
import tifffile

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.simulator.drawgaussians as dg


def clusterImages(mlist_name, title, min_size, im_max, output):
    
    i3_data = readinsight3.loadI3GoodOnly(mlist_name)

    print("Only drawing clusters with at least", min_size, "localizations.")

    rand_color = randomcolor.RandomColor()

    scale = 4
    im_size = scale * 256

    red_image = numpy.zeros((im_size, im_size))
    grn_image = numpy.zeros((im_size, im_size))
    blu_image = numpy.zeros((im_size, im_size))
    sum_image = numpy.zeros((im_size, im_size))
    
    labels = i3_data['lk']
    start = int(numpy.min(labels))
    stop = int(numpy.max(labels)) + 1
    
    num_clusters = 0
    for k in range(start,stop):
        if ((k%200)==0):
            print("Processing cluster", k)
        
        mask = (labels == k)
        csize = mask.sum()

        x = scale * i3_data['xc'][mask]
        y = scale * i3_data['yc'][mask]

        if (k == -1) or (csize < min_size):
            r = 1.0
            g = 1.0
            b = 1.0
        else:
            num_clusters += 1
            color = rand_color.generate()
            r = int(color[0][1:3], 16)/255.0
            g = int(color[0][3:5], 16)/255.0        
            b = int(color[0][5:7], 16)/255.0

        if False:
            temp = numpy.zeros((im_size, im_size))
            dg.drawGaussiansXYOnImage(temp, x, y, sigma = 1.5)

            red_image += r * temp
            grn_image += g * temp
            blu_image += b * temp
            sum_image += temp
        
        else:
            for i in range(x.size):
                yi = int(x[i])
                xi = int(y[i])
                if (xi >= 0) and (xi < im_size) and (yi >= 0) and (yi < im_size):
                    red_image[xi,yi] += r
                    grn_image[xi,yi] += g
                    blu_image[xi,yi] += b
                    sum_image[xi,yi] += 1
        
    # Some wacky normalization scheme..
    mask = (sum_image > im_max)

    red_image[mask] = (red_image[mask]/sum_image[mask]) * im_max
    grn_image[mask] = (grn_image[mask]/sum_image[mask]) * im_max
    blu_image[mask] = (blu_image[mask]/sum_image[mask]) * im_max
    sum_image[mask] = im_max

    rgb_image = numpy.zeros((im_size, im_size, 3))
    rgb_image[:,:,0] = red_image
    rgb_image[:,:,1] = grn_image
    rgb_image[:,:,2] = blu_image

    rgb_image = (255.0/im_max)*rgb_image
    sum_image = (255.0/im_max)*sum_image

    rgb_image = rgb_image.astype(numpy.uint8)
    sum_image = sum_image.astype(numpy.uint8)

    title += " (" + str(num_clusters) + ")"

    # Save RGB image.
    img1 = Image.fromarray(rgb_image, "RGB")
    try:
        draw1 = ImageDraw.Draw(img1)
        font1 = ImageFont.truetype("FreeMono.ttf", 24)
        draw1.text((2,2), title, (255,255,255), font = font1)
    except IOError:
        print("Text drawing disabled, true type font file may be missing?")
    
    img1.save(output + "_01.png")

    # Save grayscale image.
    img2 = Image.fromarray(sum_image, "L")
    try:
        draw2 = ImageDraw.Draw(img2)
        font2 = ImageFont.truetype("FreeMono.ttf", 24)
        draw2.text((2,2), title, (255,255,255), font = font2)
    except IOError:
        print("Text drawing disabled, true type font file may be missing?")
    
    img2.save(output + "_02.png")


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Make images of the identified clusters in a cell')

    parser.add_argument('--bin', dest='mlist', type=str, required=True)
    parser.add_argument('--title', dest='title', type=str, required=True)
    parser.add_argument('--min_size', dest='min_size', type=int, required=True)
    parser.add_argument('--im_max', dest='im_max', type=int, required=True)
    parser.add_argument('--image', dest='image', type=str, required=True)
    
    args = parser.parse_args()
    
    clusterImages(args.mlist, args.title, args.min_size, args.im_max, args.image)

