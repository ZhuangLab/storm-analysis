#!/usr/bin/env python
"""
Make data for testing PSF FFT.

The default tests are pretty easy as they are just relatively bright
localizations on a grid.

Hazen 10/17
"""
import storm_analysis.diagnostics.make_data_common as makeDataCommon

import storm_analysis.diagnostics.psf_fft.settings as settings


def makeData(dither = False):
    """
    Ideal camera movies.
    """
    makeDataCommon.makeDataPupilFn(settings, dither)


def makeDataRQE(cal_file = "calib.npy", dither = False):
    """
    RQE testing camera movies.
    """
    makeDataCommon.makeDataPupilFnCMOS(settings, cal_file, dither)

    
if (__name__ == "__main__"):
    makeData()
    
