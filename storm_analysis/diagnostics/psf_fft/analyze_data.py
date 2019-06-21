#!/usr/bin/env python
"""
Analyze test data using PSF FFT.

Hazen 10/17
"""
import storm_analysis.diagnostics.analyze_data_common as analyzeDataCommon

import storm_analysis.psf_fft.psffft_analysis as psffftAna


def analyzeData():
    analyzeDataCommon.analyzeData(psffftAna, "psf_fft.xml")


if (__name__ == "__main__"):
    analyzeData()
    
