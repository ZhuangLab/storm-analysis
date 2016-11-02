#!/usr/bin/env python

import numpy

import storm_analysis


def test_dbscan1():

    from storm_analysis.dbscan.dbscan_c import dbscan

    x = numpy.random.random(20)
    y = numpy.random.random(20)
    z = numpy.zeros(20)
    c = numpy.zeros(20)

    clusters = dbscan(x, y, z, c, 1.0, 10)

    for elt in clusters:
        assert (elt == 2)


if (__name__ == "__main__"):
    test_dbscan1()
    
    
