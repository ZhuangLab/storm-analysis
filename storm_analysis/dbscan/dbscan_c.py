#!/usr/bin/python
#
# 11/11
#
# Python interface to dbscan.so library.
#
# Hazen
#
# 

import ctypes
import math
import numpy
from numpy.ctypeslib import ndpointer
import os
import sys

import storm_analysis.sa_library.loadclib as loadclib

lib_dbscan = loadclib.loadCLibrary("storm_analysis.dbscan", "_dbscan")

lib_dbscan.dbscan.argtypes = [ndpointer(dtype=numpy.float32),
                              ndpointer(dtype=numpy.float32),
                              ndpointer(dtype=numpy.float32),
                              ndpointer(dtype=numpy.int32),
                              ndpointer(dtype=numpy.int32),
                              ctypes.c_int,
                              ctypes.c_float,
                              ctypes.c_int,
                              ctypes.c_int]

lib_dbscan.locClSize.argtypes = [ndpointer(dtype=numpy.int32),
                                 ndpointer(dtype=numpy.int32),
                                 ctypes.c_int,
                                 ctypes.c_int]

lib_dbscan.recategorize.argtypes = [ndpointer(dtype=numpy.int32),
                                    ndpointer(dtype=numpy.int32),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int]

#
# dbscan function.
#
# z_factor adjusts for the z resolution being about 1/2
# that of the x-y resolution.
#
# FIXME: This might be even faster (when using the kd-tree
#   approach) if the data were shuffled?
#
def dbscan(x, y, z, c, eps, min_points, z_factor = 0.5, verbose = True):
    n_peaks = x.size

    l = numpy.zeros(n_peaks, dtype = numpy.int32)

    c_x = numpy.ascontiguousarray(x.astype(numpy.float32))
    c_y = numpy.ascontiguousarray(y.astype(numpy.float32))
    c_z = numpy.ascontiguousarray(z.astype(numpy.float32))*z_factor
    c_c = numpy.ascontiguousarray(c.astype(numpy.int32))
    c_l = numpy.ascontiguousarray(l)
    lib_dbscan.dbscan(c_x,
                      c_y,
                      c_z,
                      c_c,
                      c_l,
                      n_peaks,
                      eps,
                      min_points,
                      int(verbose))

    # Print number of clusters
    if verbose:
        n_clusters_ = len(set(c_l)) - (1 if -1 in c_l else 0)
        print('Estimated number of clusters: %d' % n_clusters_)

    return c_l

#
# localizationClusterSize function.
#
# This returns the size of the cluster associated 
# with each localization.
#
def localizationClusterSize(k):
    n_peaks = k.size
    max_id = int(numpy.max(k))

    c_k = numpy.ascontiguousarray(k.astype(numpy.int32))
    c_sz = numpy.ascontiguousarray(numpy.zeros(n_peaks, dtype = numpy.int32))
    lib_dbscan.locClSize(c_sz,
                         c_k,
                         n_peaks,
                         max_id)

    return c_sz

#
# recategorize function.
#
# Note that this assumes that cluster numbers are assigned
# as by the dbscan algorithm, i.e. "good" cluster numbers
# start at 2.
#
def recategorize(k, c, min_cnts):
    n_peaks = k.size
    max_id = int(numpy.max(k))

    c_k = numpy.ascontiguousarray(k.astype(numpy.int32))
    c_c = numpy.ascontiguousarray(c.astype(numpy.int32).copy())
    lib_dbscan.recategorize(c_k,
                            c_c,
                            n_peaks,
                            max_id,
                            min_cnts)

    return c_c


