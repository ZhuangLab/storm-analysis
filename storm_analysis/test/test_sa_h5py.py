#!/usr/bin/env python
"""
Tests of our HDF5 reader, or perhaps more accurately test of our 
understanding of how to use the h5py module.
"""
import os
import storm_analysis

import storm_analysis.sa_library.sa_h5py as saH5py

def test_sa_h5py_1():
    """
    Test metadata round trip.
    """
    metadata = "<xml><field1><data1>data</data1></field></xml>"

    filename = "test_sa_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write metadata.
    with saH5py.SAH5Py(h5_name) as h5:
        h5.addMetadata(metadata)

    # Read metadata.
    with saH5py.SAH5Py(h5_name) as h5:
        print(h5.getMetadata(), type(h5.getMetadata()))

    
if (__name__ == "__main__"):
    test_sa_h5py_1()
    

    
