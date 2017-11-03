#!/usr/bin/env python
"""
Tests some basic aspects of the multifit.c / dao_fit.c libraries.
"""
import numpy

import storm_analysis.sa_library.dao_fit_c as daoFitC


def test_mfit_1():

    image = numpy.ones((40,40))
    
    mfit = daoFitC.MultiFitter()
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Check correct initialization.
    assert (mfit.getNFit() == 0)
    assert (mfit.getNFitMax() == 0)

    # Create peaks in the center of the image.
    n_peaks = 1020
    peaks = numpy.zeros((n_peaks,4))
    peaks[:,0] += 20.0
    peaks[:,1] += 20.0
    peaks[:,3] = 1.0

    # Add peaks & check size.
    mfit.newPeaks(peaks, "finder")
    assert (mfit.getNFit() == n_peaks)
    assert (mfit.getNFitMax() == 1500)

    # Add again & check size.
    peaks[:,3] = 2.0
    mfit.newPeaks(peaks, "finder")
    assert (mfit.getNFit() == 2*n_peaks)
    assert (mfit.getNFitMax() == 2500)
    
    mfit.cleanup(verbose = False)

    
if (__name__ == "__main__"):
    test_mfit_1()
    
