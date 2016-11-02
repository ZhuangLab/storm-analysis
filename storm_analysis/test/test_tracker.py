#!/usr/bin/env python

import storm_analysis


def test_tracker():
    
    # Test tracking.
    import shutil

    settings = storm_analysis.getData("test/data/test_drift.xml")
    alist_name = storm_analysis.getPathOutputTest("test_drift_alist.bin")

    # Copy mlist so that it is in the same directory as alist.
    mlist_data = storm_analysis.getData("test/data/test_drift_mlist.bin")
    mlist_output = storm_analysis.getPathOutputTest("test_drift_mlist.bin")
    shutil.copyfile(mlist_data, mlist_output)

    from storm_analysis.sa_utilities.track_average_correct import trackAverageCorrect
    trackAverageCorrect(mlist_output, alist_name, settings)


if (__name__ == "__main__"):
    test_tracker()
    
