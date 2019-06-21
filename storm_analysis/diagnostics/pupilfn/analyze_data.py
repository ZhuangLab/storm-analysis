#!/usr/bin/env python
"""
Analyze test data using Pupilfn.

Hazen 10/17
"""
import storm_analysis.diagnostics.analyze_data_common as analyzeDataCommon

import storm_analysis.pupilfn.pupilfn_analysis as pfAna


def analyzeData():
    analyzeDataCommon.analyzeData(pfAna, "pupilfn.xml")


if (__name__ == "__main__"):
    analyzeData()
    
