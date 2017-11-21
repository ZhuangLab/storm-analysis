#!/usr/bin/env python
"""
Tests for multi_plane.mp_utilities
"""

import numpy
import pickle

import storm_analysis

import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.multi_plane.mp_utilities as mpUtil

def test_load_mappings_1():
    map_test_file = storm_analysis.getPathOutputTest("map.map")

    max_ch = 4
    mappings = {}
    for i in range(1,max_ch):
        j = i
        mappings[str(i) + "_0_x"] = numpy.arange(j,j+2.5,1.0)
        j += 0.1
        mappings[str(i) + "_0_y"] = numpy.arange(j,j+2.5,1.0)
        j += 0.1
        mappings["0_" + str(i) + "_x"] = numpy.arange(j,j+2.5,1.0)
        j += 0.1
        mappings["0_" + str(i) + "_y"] = numpy.arange(j,j+2.5,1.0)

    max_ch -= 1

    with open(map_test_file, 'wb') as fp:
        pickle.dump(mappings, fp)

    mappings = {}
    [xt_0toN, yt_0toN, xt_Nto0, yt_Nto0] = mpUtil.loadMappings(map_test_file, 0)
    assert(xt_0toN[0,0] == 0.0)
    assert(yt_0toN[0,0] == 0.0)
    assert(xt_Nto0[0,0] == 0.0)
    assert(yt_Nto0[0,0] == 0.0)

    assert(abs(xt_0toN[max_ch,2]-5.2) < 1.0e-6)
    assert(abs(yt_0toN[max_ch,2]-5.3) < 1.0e-6)
    assert(abs(xt_Nto0[max_ch,2]-5.0) < 1.0e-6)
    assert(abs(yt_Nto0[max_ch,2]-5.1) < 1.0e-6)

def test_load_mappings_2():
    [xt_0toN, yt_0toN, xt_Nto0, yt_Nto0] = mpUtil.loadMappings(None, 0)
    assert(xt_0toN.shape[0] == 1)
    assert(abs(xt_0toN[0,0]) < 1.0e-6)
    assert(abs(xt_0toN[0,1] - 1.0) < 1.0e-6)
    assert(abs(xt_0toN[0,2]) < 1.0e-6)

    
if (__name__ == "__main__"):
    numpy.set_printoptions(precision=1, suppress=True)
    test_load_mappings_1()
    test_load_mappings_2()
