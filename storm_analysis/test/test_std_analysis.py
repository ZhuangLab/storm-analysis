#!/usr/bin/env python
"""
Tests of sa_utilities.std_analysis
"""
import numpy

import storm_analysis

import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_utilities.std_analysis as stdAnalysis

def test_std_analysis_1():
    """
    Test zCheck.
    """
    # Load 3D parameters.
    settings = storm_analysis.getData("test/data/test_3d_3d.xml")
    parameters = params.ParametersDAO().initFromFile(settings)

    [min_z, max_z] = parameters.getZRange()
    assert(abs(min_z + 0.5) < 1.0e-6)
    assert(abs(max_z - 0.5) < 1.0e-6)
        
    # Create HDF5 file with localizations and tracks.
    zvals = numpy.arange(-1.0, 1.05, 0.2)
    
    peaks = {"category" : numpy.ones(zvals.size, dtype = numpy.int32),
             "x" : numpy.zeros(zvals.size),
             "z" : zvals}

    h5_name = storm_analysis.getPathOutputTest("test_sa_hdf5.hdf5")
    storm_analysis.removeFile(h5_name)
    
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.setMovieInformation(256, 256, 10, "XYZZY")
        h5.addLocalizations(peaks, 1)
        h5.addTracks(peaks)

    # Run z check on the file.
    stdAnalysis.zCheck(h5_name, parameters)

    # Check track and localization categories.
    category = numpy.ones(zvals.size, dtype = numpy.int32)
    z_mask = (zvals < min_z) | (zvals > max_z)
    category[z_mask] = 9
    
    with saH5Py.SAH5Py(h5_name) as h5:
        for fnum, locs in h5.localizationsIterator(fields = ["category"]):
            assert(numpy.allclose(locs["category"], category))

        for tracks in h5.tracksIterator(fields = ["category"]):
            assert(numpy.allclose(tracks["category"], category))
    

if (__name__ == "__main__"):
    test_std_analysis_1()

    
