#!/usr/bin/env python
"""
Analyze test data using FISTA Spliner.

Hazen 11/17
"""
import storm_analysis.diagnostics.analyze_data_common as analyzeDataCommon

import storm_analysis.spliner.spline_analysis as spAna


def analyzeData():
    analyzeDataCommon.analyzeData(spAna, "fdecon.xml")
    

if (__name__ == "__main__"):
    analyzeData()
    
