#!/usr/bin/env python

import storm_analysis


def _test_frc():
    mlist_name = storm_analysis.getData("test/data/test_drift_mlist.bin")
    results_name = storm_analysis.getPathOutputTest("test_drift_frc.txt")
    
    from storm_analysis.frc.frc_calc2d import frcCalc2d
    frcCalc2d(mlist_name, results_name, False)

    
if (__name__ == "__main__"):
    _test_frc()
