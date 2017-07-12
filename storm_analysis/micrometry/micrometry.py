#!/usr/bin/env python
"""
Given two lists of localizations, returns the transform
between them.

Hazen 07/17
"""
import numpy

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3

import storm_analysis.micrometry.quads as quads


#def micrometry(quads1, 

if (__name__ == "__main__"):


    import argparse

    parser = argparse.ArgumentParser(description = 'Micrometry - ...')

    parser.add_argument('--locs1', dest='locs1', type=str, required=True,
                        help = "The name of the 'reference' localizations file")

    parser.add_argument('--locs2', dest='locs2', type=str, required=True,
                        help = "The name of the 'other' localizations file")

    args = parser.parse_args()

    print("Making quads for the 'reference' data.")
    i3_data1 = readinsight3.loadI3File(args.locs1)

    i3_data1 = i3dtype.maskData(i3_data1, (i3_data1['fr'] == 1))
    x1 = i3_data1['xc']
    y1 = i3_data1['yc']
    h1 = i3_data1['h']


    quads1 = quads.makeQuads(x1, y1, h1, min_size = 5.0, max_size = 100.0)

    print("Making quads for the 'other' data.")
    i3_data2 = readinsight3.loadI3File(args.locs2)

    i3_data2 = i3dtype.maskData(i3_data2, (i3_data2['fr'] == 1))
    x2 = i3_data2['xc']
    y2 = i3_data2['yc']
    h2 = i3_data2['h']

    quads2 = quads.makeQuads(x2, y2, h2, min_size = 5.0, max_size = 100.0)    
    print(len(quads1), len(quads2))

    print("Comparing.")
    matches = 0
    for q1 in quads1:
        for q2 in quads2:
            if q1.isMatch(q2):
                print(q1)
                print(q2)
                print("Transform:", q1.getTransform(q2))
                print()
                matches += 1

    print("Found", matches)
