#!/usr/bin/env python
"""
Tests some basic aspects of the multifit.c / dao_fit.c libraries.
"""

import numpy
import tifffile

import storm_analysis.sa_library.dao_fit_c as daoFitC
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.simulator.draw_gaussians_c as dg


def test_mfit_1():
    """
    Test initialization & growing peak array size.
    """
    image = numpy.ones((40,40))
    
    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Check correct initialization.
    assert (mfit.getNFit() == 0)
    assert (mfit.getNFitMax() == 0)

    # Create peaks in the center of the image.
    n_peaks = 1020
    peaks = {"x" : numpy.ones(n_peaks) * 20.0,
             "y" : numpy.ones(n_peaks) * 20.0,
             "z" : numpy.zeros(n_peaks),
             "sigma" : numpy.ones(n_peaks)}

    # Add peaks & check size.
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == n_peaks)
    assert (mfit.getNFitMax() == 1500)

    # Add again & check size.
    peaks["sigma"][:] = 2.0
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == 2*n_peaks)
    assert (mfit.getNFitMax() == 2500)

    # Check some peak values.
    w = mfit.getPeakProperty("xsigma")
    assert (abs(w[0]-1.0) < 1.0e-6)
    assert (abs(w[n_peaks]-2.0) < 1.0e-6)
    
    mfit.cleanup(verbose = False)

def test_mfit_2():
    """
    Test initialization & growing with error peak removal.
    """
    image = numpy.ones((40,40))

    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Create ERROR peaks (too close to edge).
    n_peaks = 5
    peaks = {"x" : numpy.zeros(n_peaks),
             "y" : numpy.zeros(n_peaks),
             "z" : numpy.zeros(n_peaks),
             "sigma" : numpy.ones(n_peaks)}

    # Add peaks & check size.
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == 5)
    assert (mfit.getNFitMax() == 500)

    # Add good peaks. ERROR peaks should be removed.
    n_peaks = 500
    peaks = {"x" : numpy.ones(n_peaks) * 20.0,
             "y" : numpy.ones(n_peaks) * 20.0,
             "z" : numpy.zeros(n_peaks),
             "sigma" : numpy.ones(n_peaks) * 2.0}
    mfit.newPeaks(peaks, "testing")

    assert (mfit.getNFit() == n_peaks)
    assert (mfit.getNFitMax() == 1000)

    # Check for correct peak values.
    w = mfit.getPeakProperty("xsigma")
    assert (abs(w[0]-2.0) < 1.0e-6)

    mfit.cleanup(verbose = False)

def test_mfit_3():
    """
    Test error peak removal.
    """
    image = numpy.ones((40,40))
    
    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Add good peaks.
    n_peaks = 2
    peaks = {"x" : numpy.ones(n_peaks) * 20.0,
             "y" : numpy.ones(n_peaks) * 20.0,
             "z" : numpy.zeros(n_peaks),
             "sigma" : numpy.ones(n_peaks)}
    mfit.newPeaks(peaks, "testing")

    # Check that no peaks are removed.
    mfit.removeErrorPeaks()
    assert (mfit.getNFit() == 2)

    # Add error peaks.
    n_peaks = 2
    peaks = {"x" : numpy.zeros(n_peaks),
             "y" : numpy.zeros(n_peaks),
             "z" : numpy.zeros(n_peaks),
             "sigma" : numpy.ones(n_peaks) * 2.0}
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == 4)

    # Add more good peaks.
    n_peaks = 2
    peaks = {"x" : numpy.ones(n_peaks) * 20.0,
             "y" : numpy.ones(n_peaks) * 20.0,
             "z" : numpy.zeros(n_peaks),
             "sigma" : numpy.ones(n_peaks)}
    mfit.newPeaks(peaks, "testing")
    assert (mfit.getNFit() == 6)

    # Check that two peaks are removed.
    mfit.removeErrorPeaks()
    assert (mfit.getNFit() == 4)

    # Check that the right peaks were removed.
    w = mfit.getPeakProperty("xsigma")
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
    
    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(background)

    peaks = {"x" : numpy.array([50.0, 54.0]),
             "y" : numpy.array([50.0, 50.0]),
             "z" : numpy.array([0.0, 0.0]),
             "sigma" : numpy.array([sigma, sigma])}
    
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
    
    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(background)

    peaks = {"x" : numpy.array([50.0, 54.0]),
             "y" : numpy.array([50.0, 50.0]),
             "z" : numpy.array([0.0, 0.0]),
             "sigma" : numpy.array([sigma, sigma])}

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
        assert (abs(x[i] - peaks["x"][i]) < 1.0e-2)
        assert (abs(y[i] - peaks["y"][i]) < 1.0e-2)
    
    # Check peak w.
    w = mfit.getPeakProperty("xsigma")
    for i in range(w.size):
        assert (abs((w[i] - sigma)/sigma) < 0.02)
    
    mfit.cleanup(verbose = False)
    
def test_mfit_6():
    """
    Test marking peak status.
    """
    image = numpy.ones((40,40))
    
    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Add good peaks.
    n_peaks = 5
    peaks = {"x" : numpy.ones(n_peaks) * 20.0,
             "y" : numpy.ones(n_peaks) * 20.0,
             "z" : numpy.zeros(n_peaks),
             "sigma" : 0.1*numpy.arange(n_peaks) + 1.0}
    mfit.newPeaks(peaks, "testing")

    # Check that no peaks are removed.
    mfit.removeErrorPeaks()
    assert (mfit.getNFit() == 5)

    # Mark every other peak as bad.
    status = mfit.getPeakProperty("status")
    status[1] = iaUtilsC.ERROR
    status[3] = iaUtilsC.ERROR
    mfit.setPeakStatus(status)
    
    # Check that two peaks were removed.
    mfit.removeErrorPeaks()
    assert (mfit.getNFit() == 3)

    # Check that the right peaks were removed.
    w = mfit.getPeakProperty("xsigma")
    assert numpy.allclose(w, numpy.array([1.0, 1.2, 1.4]))

    mfit.cleanup(verbose = False)

def test_mfit_7():
    """
    Test removing RUNNING peaks.
    """
    image = numpy.ones((40,40))
    
    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Add good peaks.
    n_peaks = 2
    peaks = {"x" : numpy.ones(n_peaks) * 20.0,
             "y" : numpy.ones(n_peaks) * 20.0,
             "z" : numpy.zeros(n_peaks),
             "sigma" : numpy.ones(n_peaks) * 1.0}
    peaks["sigma"][1] = 2.0
    mfit.newPeaks(peaks, "testing")

    # Mark peaks as converged.
    status = mfit.getPeakProperty("status")
    status[:] = iaUtilsC.CONVERGED
    mfit.setPeakStatus(status)
    
    # Check that no peaks are removed.
    mfit.removeRunningPeaks()
    assert (mfit.getNFit() == 2)

    # Mark the first peak as running.
    status = mfit.getPeakProperty("status")
    status[0] = iaUtilsC.RUNNING
    mfit.setPeakStatus(status)
    
    # Check that one peak was removed.
    mfit.removeRunningPeaks()
    assert (mfit.getNFit() == 1)

    # Check that the right peak was removed.
    w = mfit.getPeakProperty("xsigma")
    for i in range(w.size):
        assert (abs(w[i] - 2.0) < 1.0e-6)

    mfit.cleanup(verbose = False)    

def test_mfit_8():
    """
    Test 'pre-specified' peak locations addition (text).
    """
    image = numpy.ones((40,40))
    
    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Add peaks.
    peaks = {"background" : numpy.array([10.0, 20.0]),
             "height" : numpy.array([11.0, 21.0]),
             "x" : numpy.array([12.0, 22.0]),
             "xsigma" : numpy.array([1.0, 2.0]),
             "y" : numpy.array([14.0, 24.0]),
             "ysigma" : numpy.array([3.0, 4.0]),
             "z" : numpy.array([16.0, 26.0])}
    mfit.newPeaks(peaks, "text")

    # Round trip verification.
    for pname in peaks:
        pvals = peaks[pname]
        mvals = mfit.getPeakProperty(pname)
        for i in range(pvals.size):
            assert(abs(pvals[i] - mvals[i]) < 1.0e-6)

    mfit.cleanup(verbose = False)

def test_mfit_9():
    """
    Test peak significance calculation.
    """
    height = 10.0
    sigma = 1.5
    x_size = 100
    y_size = 120
    background = numpy.zeros((x_size, y_size)) + 10.0
    image = dg.drawGaussians((x_size, y_size),
                             numpy.array([[50.0, 30.0, height, sigma, sigma],
                                          [50.0, 50.0, 2.0*height, sigma, sigma],
                                          [50.0, 70.0, 3.0*height, sigma, sigma]]))
    image += background

    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(background)

    peaks = {"x" : numpy.array([30.0, 50.0, 70.0]),
             "y" : numpy.array([50.0, 50.0, 50.0]),
             "z" : numpy.array([0.0, 0.0, 0.0]),
             "sigma" : numpy.array([sigma, sigma, sigma])}

    mfit.newPeaks(peaks, "finder")

    if False:
        with tifffile.TiffWriter("test_mfit_9.tif") as tf:
            tf.save(image.astype(numpy.float32))
    
    sig = mfit.getPeakProperty("significance")
    assert(abs((sig[1] - sig[0]) - (sig[2] - sig[1])) < 0.1)

def test_mfit_10():
    """
    Test 'pre-specified' peak locations addition (hdf5).
    """
    image = numpy.ones((40,40))
    
    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(image)

    # Add peaks.
    peaks = {"background" : numpy.array([10.0, 20.0]),
             "height" : numpy.array([11.0, 21.0]),
             "x" : numpy.array([12.0, 22.0]),
             "xsigma" : numpy.array([1.0, 2.0]),
             "y" : numpy.array([14.0, 24.0]),
             "ysigma" : numpy.array([3.0, 4.0]),
             "z" : numpy.array([16.0, 26.0])}
    mfit.newPeaks(peaks, "hdf5")

    # Round trip verification.
    for pname in peaks:
        pvals = peaks[pname]
        mvals = mfit.getPeakProperty(pname)
        for i in range(pvals.size):
            assert(abs(pvals[i] - mvals[i]) < 1.0e-6)

    mfit.cleanup(verbose = False)
    
    
if (__name__ == "__main__"):
    test_mfit_1()
    test_mfit_2()
    test_mfit_3()
    test_mfit_4()
    test_mfit_5()
    test_mfit_6()
    test_mfit_7()
    test_mfit_8()
    test_mfit_9()
    test_mfit_10()
    
