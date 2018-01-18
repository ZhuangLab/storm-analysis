#!/usr/bin/env python
import shutil

import storm_analysis
import storm_analysis.sa_library.parameters as params
import storm_analysis.rcc.rcc_drift_correction as rccDriftCorrection

import storm_analysis.test.verifications as veri
    
def _test_rcc_1():
    """
    Test RCC drift correction.
    """
    # Calculate drift correction.
    param_name = storm_analysis.getData("test/data/test_drift.xml")    
    parameters = params.ParametersCommon().initFromFile(param_name)
    
    data_name = storm_analysis.getData("test/data/test_drift.hdf5")
    h5_name = storm_analysis.getPathOutputTest("test_rcc_hdf5.hdf5")
    
    # Make a copy of the original as it will get modified and then
    # git will pick this up.
    shutil.copyfile(data_name, h5_name)
    
    drift_output = storm_analysis.getPathOutputTest("test_drift_drift.txt")

    [min_z, max_z] = parameters.getZRange()
    rccDriftCorrection.rccDriftCorrection(h5_name,
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

    # These thresholds are somewhat arbitrary.
    if (diffs[1] > 0.1) or (diffs[2] > 0.1):
        raise Exception("XY drift correction error.")

    if (diffs[3] > 0.03):
        raise Exception("Z drift correction error.")    

if (__name__ == "__main__"):
    _test_rcc_1()
