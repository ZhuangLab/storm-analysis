#!/usr/bin/env python

import storm_analysis

import storm_analysis.test.verify_drift_correction as vDC

def test_rcc():

    # Calculate drift correction.
    mlist_name = storm_analysis.getData("test/data/test_drift_mlist.bin")
    drift_name = storm_analysis.getPathOutputTest("test_rcc_drift.txt")

    from storm_analysis.rcc.rcc_drift_correction import rccDriftCorrection
    rccDriftCorrection(mlist_name, drift_name, 1000, 1, True, False)

    # Verify results.
    diffs = vDC.verifyDriftCorrection(storm_analysis.getData("test/data/test_drift.txt"), drift_name)
    
    if (diffs[0] > 0.1):
        raise Exception("Frame numbers do not match.")

    # These thresholds are somewhat arbitrary.
    if (diffs[1] > 0.1) or (diffs[2] > 0.1):
        raise Exception("XY drift correction error.")

    if (diffs[3] > 30.0):
        raise Exception("Z drift correction error.")    

if (__name__ == "__main__"):
    test_rcc()
