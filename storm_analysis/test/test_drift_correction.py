#!/usr/bin/env python
import numpy
import os

import storm_analysis
import storm_analysis.sa_library.drift_utilities as driftUtils
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_utilities.xyz_drift_correction as xyzDriftCorrection

import storm_analysis.test.verifications as veri

def test_drift_correction_1():
    """
    This tests the whole process.
    """
    # Calculate drift correction.
    param_name = storm_analysis.getData("test/data/test_drift.xml")    
    parameters = params.ParametersCommon().initFromFile(param_name)

    bin_name = storm_analysis.getData("test/data/test_drift.hdf5")
    drift_output = storm_analysis.getPathOutputTest("test_drift_drift.txt")

    [min_z, max_z] = parameters.getZRange()
    xyzDriftCorrection.xyzDriftCorrection(bin_name,
                                          drift_output,
                                          parameters.getAttr("frame_step"),
                                          parameters.getAttr("d_scale"),
                                          min_z,
                                          max_z,
                                          True)

    # Verify results.
    diffs = veri.verifyDriftCorrection(storm_analysis.getData("test/data/test_drift.txt"),
                                       drift_output)
    
    if (diffs[0] > 0.1):
        raise Exception("Frame numbers do not match.")

    # These thresholds are somewhat arbitrary, 0.1 pixel maximum error in X/Y, 30nm in Z.
    if (diffs[1] > 0.1) or (diffs[2] > 0.1):
        raise Exception("XY drift correction error.")

    if (diffs[3] > 0.03):
        raise Exception("Z drift correction error.")

def test_drift_correction_2():
    """
    Test interpolation.
    """
    x_vals = numpy.array([3,6,9,12])
    y_vals = numpy.array([1,2,2,1])
    int_y = driftUtils.interpolateData(x_vals, y_vals, 16)

    exp_y = numpy.array([0.0, 1.0/3.0, 2.0/3.0, 1.0, 4.0/3.0, 5.0/3.0, 2.0, 2.0,
                         2.0, 2.0, 5.0/3.0, 4.0/3.0, 1.0, 2.0/3.0, 1.0/3.0, 0.0])

    assert(numpy.allclose(int_y, exp_y))

    if False:
        import matplotlib
        import matplotlib.pyplot as pyplot

        int_x = numpy.arange(16)
        pyplot.plot(int_x, int_y)
        pyplot.show()

def test_drift_correction_3():
    """
    Test handling of files with no localizations.
    """
    filename = "test_dc_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)
    
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.setMovieProperties(128, 128, 10000, "XYZZY")

    drift_output = storm_analysis.getPathOutputTest("test_drift_drift.txt")
    
    xyzDriftCorrection.xyzDriftCorrection(h5_name,
                                          drift_output,
                                          500,
                                          2,
                                          -0.5,
                                          0.5,
                                          False)

    drift_data = numpy.loadtxt(drift_output)
    assert(numpy.allclose(drift_data[:,1], numpy.zeros(drift_data.shape[0])))

def test_drift_correction_4():
    """
    Test handling of very short files.
    """
    peaks = {"x" : numpy.zeros(10),
             "y" : numpy.ones(10)}

    filename = "test_dc_hdf5.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)
    
    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.setMovieProperties(128, 128, 100, "XYZZY")
        h5.addLocalizations(peaks, 0)
        h5.addLocalizations(peaks, 2)

    drift_output = storm_analysis.getPathOutputTest("test_drift_drift.txt")
    
    xyzDriftCorrection.xyzDriftCorrection(h5_name,
                                          drift_output,
                                          500,
                                          2,
                                          -0.5,
                                          0.5,
                                          False)

    drift_data = numpy.loadtxt(drift_output)
    assert(numpy.allclose(drift_data[:,1], numpy.zeros(drift_data.shape[0])))
    
if (__name__ == "__main__"):
    test_drift_correction_1()
    test_drift_correction_2()
    test_drift_correction_3()
    test_drift_correction_4()
