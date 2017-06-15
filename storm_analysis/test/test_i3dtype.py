#!/usr/bin/env python
"""
Tests for sa_library.i3dtype
"""
import numpy

import storm_analysis.sa_library.i3dtype as i3dtype


def test_i3dtype_1():
    x_size = 100
    y_size = 100
    frame = 10
    nm_per_pixel = 100.0
    
    data_in = i3dtype.createDefaultI3Data(10)
    i3dtype.posSet(data_in, 'x', numpy.arange(10))
    i3dtype.posSet(data_in, 'y', numpy.arange(10) + 30.0)
    i3dtype.posSet(data_in, 'z', numpy.arange(10) + 60.0)
    i3dtype.setI3Field(data_in, 'fr', frame)

    peaks = i3dtype.convertToMultiFit(data_in, x_size, y_size, frame, nm_per_pixel)
    data_out = i3dtype.createFromMultiFit(peaks, x_size, y_size, frame, nm_per_pixel)

    fields = ['x', 'ax', 'w']
    for i in range(10):
        for field in fields:
            assert(abs(data_in[field][i] - data_out[field][i]) < 1.0e-6)
    
def test_i3dtype_2():
    i3_locs = i3dtype.createDefaultI3Data(10)
    i3dtype.posSet(i3_locs, 'x', numpy.arange(10))

    for i in range(10):
        assert(abs(i3_locs["x"][i] - i) < 1.0e-6)
        assert(abs(i3_locs["xc"][i] - i) < 1.0e-6)

def test_i3dtype_3():
    i3_locs = i3dtype.createDefaultI3Data(10)
    i3dtype.posSet(i3_locs, 'x', 10.0)

    for i in range(10):
        assert(abs(i3_locs["x"][i] - 10.0) < 1.0e-6)
        assert(abs(i3_locs["xc"][i] - 10.0) < 1.0e-6)
        
    
if (__name__ == "__main__"):
    test_i3dtype_1()
    test_i3dtype_2()
    test_i3dtype_3()
    
