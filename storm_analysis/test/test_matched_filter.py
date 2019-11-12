#!/usr/bin/env python
"""
Tests for sa_library.matched_filter_c.py

Hazen 06/17
"""
import numpy

import storm_analysis.sa_library.matched_filter_c as matchedFilterC
import storm_analysis.sa_library.recenter_psf as recenterPSF
import storm_analysis.simulator.draw_gaussians_c as dg


def test_matched_filter1():
    """
    Verify that the filter results are normalized correctly.
    """
    x_size = 80
    y_size = 90

    objects = numpy.zeros((1, 5))

    # Make filter with unit sum.
    objects[0,:] = [x_size/2, y_size/2, 1.0, 1.0, 1.0]
    psf = dg.drawGaussians((x_size, y_size), objects)
    psf = psf/numpy.sum(psf)
    flt = matchedFilterC.MatchedFilter(psf)

    for i in range(1,5):
        image = numpy.zeros((x_size, y_size))
        image[int(x_size/2), int(y_size/2)] = float(i)
        conv = flt.convolve(image)
        assert(abs(numpy.sum(image) - numpy.sum(conv)) < 1.0e-6)

    flt.cleanup()
    
        
def test_matched_filter2():
    """
    Test recovering the original height from the convolved image.

    FIXME: Not sure this test is relevant anymore as we no longer use
           this approach for initializing peak heights.
    """
    x_size = 80
    y_size = 90

    objects = numpy.zeros((1, 5))

    # Make filter with unit sum.    
    objects[0,:] = [x_size/2, y_size/2, 1.0, 2.0, 2.0]
    psf = dg.drawGaussians((x_size, y_size), objects)
    psf_norm = psf/numpy.sum(psf)
    flt = matchedFilterC.MatchedFilter(psf_norm)

    rescale = 1.0/numpy.sum(psf * psf_norm)

    # Create object with height 10 and the same shape as the filter.
    height = 10.0
    objects[0,:] = [x_size/2, y_size/2, height, 2.0, 2.0]
    image = dg.drawGaussians((x_size, y_size), objects)

    # Apply filter.
    conv = flt.convolve(image)

    # Verify that final height is 'close enough'.
    assert (abs(numpy.amax(conv) * rescale - height)/height < 1.0e-2)

    flt.cleanup()

    
def test_matched_filter3():
    """
    Test memoization.
    """
    x_size = 40
    y_size = 50

    objects = numpy.zeros((1, 5))

    # Make filter with unit sum.    
    objects[0,:] = [x_size/2, y_size/2, 1.0, 2.0, 2.0]
    psf = dg.drawGaussians((x_size, y_size), objects)
    psf_norm = psf/numpy.sum(psf)
    flt = matchedFilterC.MatchedFilter(psf_norm, memoize = True, max_diff = 0.1)

    # Convolve first image.
    image = numpy.zeros((x_size, y_size))
    image[int(x_size/2),int(y_size/2)] = 1.0

    res1 = flt.convolve(image)

    # Repeat with approximately the same image.
    image = numpy.zeros((x_size, y_size))
    image[int(x_size/2),int(y_size/2)] = 1.05

    res2 = flt.convolve(image)
    assert(abs(numpy.max(res2 - res1)) < 1.0e-12)

    # Repeat with a different image.
    image = numpy.zeros((x_size, y_size))
    image[int(x_size/2),int(y_size/2)] = 1.1

    res2 = flt.convolve(image)
    assert(abs(numpy.max(res2 - res1)) > 1.0e-12)

    # Repeat with original image to verify update.
    image = numpy.zeros((x_size, y_size))
    image[int(x_size/2),int(y_size/2)] = 1.0

    res1 = flt.convolve(image)
    assert(abs(numpy.max(res2 - res1)) > 1.0e-12)

    # Repeat with exactly the same image.
    image = numpy.zeros((x_size, y_size))
    image[int(x_size/2),int(y_size/2)] = 1.0

    res2 = flt.convolve(image)
    assert(abs(numpy.max(res2 - res1)) < 1.0e-12)

    flt.cleanup()

    
def test_matched_filter4():
    """
    Verify that fftw_estimate has no effect on the results.
    """
    x_size = 80
    y_size = 90

    objects = numpy.zeros((1, 5))

    # Make filter with unit sum.
    objects[0,:] = [x_size/2, y_size/2, 1.0, 1.0, 1.0]
    psf = dg.drawGaussians((x_size, y_size), objects)
    psf = psf/numpy.sum(psf)
    flt1 = matchedFilterC.MatchedFilter(psf, fftw_estimate = False)
    flt2 = matchedFilterC.MatchedFilter(psf, fftw_estimate = True)

    for i in range(1,5):
        image = numpy.zeros((x_size, y_size))
        image[int(x_size/2), int(y_size/2)] = float(i)
        conv1 = flt1.convolve(image)
        conv2 = flt2.convolve(image)
        
        assert(numpy.allclose(conv1, conv2))

    flt1.cleanup()
    flt2.cleanup()


def test_matched_filter5():
    """
    Verify that the filter results are correct.
    """
    x_size = 80
    y_size = 90

    objects = numpy.zeros((1, 5))

    # Make filter with unit sum.
    objects[0,:] = [x_size/2, y_size/2, 1.0, 1.0, 1.0]
    psf = dg.drawGaussians((x_size, y_size), objects)
    psf = psf/numpy.sum(psf)
    flt = matchedFilterC.MatchedFilter(psf)

    # Make test image.
    image = numpy.zeros((x_size, y_size))
    image[int(x_size/2), int(y_size/2)] = float(100)

    mf_conv = flt.convolve(image)

    t1 = numpy.fft.fft2(recenterPSF.recenterPSF(psf))
    t2 = numpy.fft.fft2(image)
    np_conv = numpy.real(numpy.fft.ifft2(t1*t2))

    assert(numpy.allclose(mf_conv, np_conv))

    flt.cleanup()
    

if (__name__ == "__main__"):
    test_matched_filter5()

