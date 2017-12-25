#!/usr/bin/env python

import numpy

import storm_analysis
import storm_analysis.sa_library.drift_utilities as driftUtils
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.xyz_drift_correction as xyzDriftCorrection

import storm_analysis.test.verifications as veri

def test_drift_correction_1():
    """
    This tests the whole process.
    """
    # Calculate drift correction.
    param_name = storm_analysis.getData("test/data/test_drift.xml")    
    parameters = params.ParametersCommon().initFromFile(param_name)

    mlist_name = storm_analysis.getData("test/data/test_drift.hdf5")
    drift_output = storm_analysis.getPathOutputTest("test_drift_drift.txt")

    [min_z, max_z] = parameters.getZRange()
    xyzDriftCorrection.xyzDriftCorrection(mlist_name,
                                          drift_output,
                                          parameters.getAttr("frame_step"),
                                          parameters.getAttr("d_scale"),
                                          min_z,
                                          max_z,
                                          True)

    # Verify results.
    diffs = veri.verifyDriftCorrection(storm_analysis.getData("test/data/test_drift.txt"), drift_output)
    
    if (diffs[0] > 0.1):
        raise Exception("Frame numbers do not match.")

    # These thresholds are somewhat arbitrary.
    if (diffs[1] > 0.1) or (diffs[2] > 0.1):
        raise Exception("XY drift correction error.")

    if (diffs[3] > 30.0):
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

    
if (__name__ == "__main__"):
    test_drift_correction_1()
    test_drift_correction_2()
    
