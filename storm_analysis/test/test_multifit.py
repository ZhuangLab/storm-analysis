#!/usr/bin/env python
"""
Tests some basic aspects of the multifit.c / dao_fit.c libraries.
"""
import numpy
import tifffile

import storm_analysis.sa_library.dao_fit_c as daoFitC
import storm_analysis.simulator.draw_gaussians_c as dg


def test_mfit_1():
    """
    Test initialization & growing peak array size.
    """
    image = numpy.ones((40,40))
    
    mfit = daoFitC.MultiFitter2D()
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Check correct initialization.
    assert (mfit.getNFit() == 0)
    assert (mfit.getNFitMax() == 0)

    # Create peaks in the center of the image.
    n_peaks = 1020
    peaks = numpy.zeros((n_peaks,4))
    peaks[:,0] = 20.0
    peaks[:,1] = 20.0
    peaks[:,3] = 1.0

    # Add peaks & check size.
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == n_peaks)
    assert (mfit.getNFitMax() == 1500)

    # Add again & check size.
    peaks[:,3] = 2.0
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == 2*n_peaks)
    assert (mfit.getNFitMax() == 2500)

    # Check some peak values.
    w = mfit.getPeakProperty("xwidth")
    assert (abs(w[0]-1.0) < 1.0e-6)
    assert (abs(w[n_peaks]-2.0) < 1.0e-6)
    
    mfit.cleanup(verbose = False)

def test_mfit_2():
    """
    Test initialization & growing with error peak removal.
    """
    image = numpy.ones((40,40))

    mfit = daoFitC.MultiFitter2D()
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Create ERROR peaks (too close to edge).
    peaks = numpy.zeros((5,4))
    peaks[:,3] = 1.0

    # Add peaks & check size.
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == 5)
    assert (mfit.getNFitMax() == 500)

    # Add good peaks. ERROR peaks should be removed.
    n_peaks = 500
    peaks = numpy.zeros((n_peaks,4))
    peaks[:,0] = 20.0
    peaks[:,1] = 20.0
    peaks[:,3] = 2.0
    mfit.newPeaks(peaks, "testing")

    assert (mfit.getNFit() == n_peaks)
    assert (mfit.getNFitMax() == 1000)

    # Check for correct peak values.
    w = mfit.getPeakProperty("xwidth")
    assert (abs(w[0]-2.0) < 1.0e-6)

    mfit.cleanup(verbose = False)

def test_mfit_3():
    """
    Test error peak removal.
    """
    image = numpy.ones((40,40))
    
    mfit = daoFitC.MultiFitter2D()
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Add good peaks.
    peaks = numpy.zeros((2,4))
    peaks[:,0] = 20.0
    peaks[:,1] = 20.0
    peaks[:,3] = 1.0    
    mfit.newPeaks(peaks, "testing")

    # Check that no peaks are removed.
    mfit.removeErrorPeaks()
    assert (mfit.getNFit() == 2)

    # Add error peaks.
    peaks = numpy.zeros((2,4))
    peaks[:,3] = 2.0
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == 4)

    # Add more good peaks.
    peaks = numpy.zeros((2,4))
    peaks[:,0] = 20.0
    peaks[:,1] = 20.0
    peaks[:,3] = 1.0
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == 6)

    # Check that two peaks are removed.
    mfit.removeErrorPeaks()
    assert (mfit.getNFit() == 4)

    # Check that the right peaks were removed.
    w = mfit.getPeakProperty("xwidth")
    for i in range(w.size):
        assert (abs(w[i] - 1.0) < 1.0e-6)

    mfit.cleanup(verbose = False)

def test_mfit_4():
    """
    Test height and background initialization.
    """
    height = 20.0
    sigma = 1.5
    x_size = 100
    y_size = 120
    background = numpy.zeros((x_size, y_size)) + 10.0
    image = dg.drawGaussians((x_size, y_size),
                             numpy.array([[50.0, 50.0, height, sigma, sigma],
                                          [50.0, 54.0, height, sigma, sigma]]))
    image += background
    
    mfit = daoFitC.MultiFitter2D()
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(background)

    peaks = numpy.array([[50.0, 50.0, 0.0, sigma],
                         [54.0, 50.0, 0.0, sigma]])

    mfit.newPeaks(peaks, "finder")

    if False:
        with tifffile.TiffWriter("test_mfit_4.tif") as tf:
            tf.save((image-background).astype(numpy.float32))
            tf.save(mfit.getFitImage().astype(numpy.float32))

    # Check height.
    h = mfit.getPeakProperty("height")
    for i in range(h.size):
        assert (abs(h[i] - height)/height < 0.1)

    # Check background.
    bg = mfit.getPeakProperty("background")
    for i in range(bg.size):
        assert (abs(bg[i] - 10.0) < 1.0e-6)

    mfit.cleanup(verbose = False)
        
def test_mfit_5():
    """
    Test initialization and fitting.
    """
    height = 20.0
    sigma = 1.5
    x_size = 100
    y_size = 120
    background = numpy.zeros((x_size, y_size)) + 10.0
    image = dg.drawGaussians((x_size, y_size),
                             numpy.array([[50.0, 50.0, height, sigma, sigma],
                                          [50.0, 54.0, height, sigma, sigma]]))
    image += background
    
    mfit = daoFitC.MultiFitter2D()
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(background)

    peaks = numpy.array([[50.0, 50.0, 0.0, sigma],
                         [54.0, 50.0, 0.0, sigma]])

    mfit.newPeaks(peaks, "finder")
    mfit.doFit()

    if False:
        with tifffile.TiffWriter("test_mfit_5.tif") as tf:
            tf.save((image-background).astype(numpy.float32))
            tf.save(mfit.getFitImage().astype(numpy.float32))
    
    # Check peak x,y.
    x = mfit.getPeakProperty("x")
    y = mfit.getPeakProperty("y")
    for i in range(x.size):
        assert (abs(x[i] - peaks[i,0]) < 1.0e-3)
        assert (abs(y[i] - peaks[i,1]) < 1.0e-3)
    
    # Check peak w.
    w = mfit.getPeakProperty("xwidth")
    for i in range(w.size):
        assert (abs((w[i] - sigma)/sigma) < 0.02)
    
    mfit.cleanup(verbose = False)
    

if (__name__ == "__main__"):
    test_mfit_1()
    test_mfit_2()
    test_mfit_3()
    test_mfit_4()
    test_mfit_5()

