#!/usr/bin/env python

import numpy

import storm_analysis
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.xyz_drift_correction as xyzDriftCorrection

import storm_analysis.test.verifications as veri

def test_drift_correction():

    # Calculate drift correction.
    param_name = storm_analysis.getData("test/data/test_drift.xml")    
    parameters = params.ParametersCommon().initFromFile(param_name)

    mlist_name = storm_analysis.getData("test/data/test_drift_mlist.bin")
    drift_output = storm_analysis.getPathOutputTest("test_drift_drift.txt")

    [min_z, max_z] = parameters.getZRange()
    xyzDriftCorrection.xyzDriftCorrection(mlist_name,
                                          drift_output,
                                          parameters.getAttr("frame_step"),
                                          parameters.getAttr("d_scale"),
                                          1000.0 * min_z,
                                          1000.0 * max_z,
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


if (__name__ == "__main__"):
    test_drift_correction()
    
