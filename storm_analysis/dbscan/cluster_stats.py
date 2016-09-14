#!/usr/bin/env python
#
# Prints some stats of the clusters that contain
# at least the user spec'd number of localizations.
#
# Hazen 11/11
#

import math
import numpy
import sys

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3

if(__name__ == "__main__"):
    
    # Load the data.
    pix_to_nm = 160.0
    
    i3_data_in = readinsight3.loadI3GoodOnly(sys.argv[1])
    stats_fp = open(sys.argv[1][:-8] + "stats.txt", "w")
    header = ["cluster", "cat", "size", "x-center(nm)", "y-center(nm)", "z-center(nm)", "size-x(nm)", "size-y(nm)", "size-z(nm)", "rg"]
    stats_fp.write(" ".join(header) + "\n")
    
    # Remove category zero localizations.
    i3_data = i3dtype.maskData(i3_data_in, (i3_data_in['c'] != 0))

    # Calculate cluster stats.
    labels = i3_data['lk']
    for k in set(labels):
        if (k>=2):
            mask = (labels == k)
            csize = mask.sum()
            if (csize > int(sys.argv[2])):
                x = pix_to_nm*i3_data['xc'][mask]
                y = pix_to_nm*i3_data['yc'][mask]
                z = i3_data['zc'][mask]
                c = i3_data['c'][mask]

                # Calculate size in x, y, z.
                sx = numpy.max(x) - numpy.min(x)
                sy = numpy.max(y) - numpy.min(y)
                sz = numpy.max(z) - numpy.min(z)

                # Calculate radius of gyration.
                cx = numpy.mean(x)
                cy = numpy.mean(y)

                rx = x - cx
                ry = y - cy

                if 0:
                    cz = numpy.mean(z)
                    rz = z - cz
                    rg = math.sqrt(numpy.sum(rx*rx + ry*ry + rz*rz) / float(x.size))
                else:
                    rg = math.sqrt(numpy.sum(rx*rx + ry*ry) / float(x.size))

                print("Cluster:", k, x.size, "localizations")
                stats = map(str, [k, c[0], csize, numpy.mean(x), numpy.mean(y), numpy.mean(z), sx, sy, sz, rg])
                stats_fp.write(" ".join(stats) + "\n")

    stats_fp.close()
