#!/usr/bin/env python
"""
Tests some aspect of ia_utilities.
"""
import numpy
import tifffile

import storm_analysis.sa_library.dao_fit_c as daoFitC
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.simulator.draw_gaussians_c as dg


def imagesCopy(images):
    images_copy = []
    for image in images:
        images_copy.append(numpy.copy(image))
    return images_copy

        
def test_ia_util_1():
    """
    Test finding peaks in an empty image.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(images)
    assert(x.size == 0)

def test_ia_util_2():
    """
    Test finding peaks in an image.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    # Above threshold, unequal height.
    images[0][10,11] = 1.1
    images[0][10,12] = 1.5

    # Above threshold, equal height.
    images[0][15,8] = 1.5
    images[0][15,9] = 1.5

    # Below threshold.
    images[0][20,11] = 0.5
    
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z, h] = mxf.findMaxima(images, want_height = True)

    assert(x.size == 2)
    for i in range(z.size):
        assert (abs(z[i] - z_values[0]) < 1.0e-6)
        assert (abs(h[i] - 1.5) < 1.0e-6)

def test_ia_util_3():
    """
    Test agreement with fitting regarding orientation.
    """
    height = 20.0
    sigma = 1.5
    x_size = 100
    y_size = 120
    background = numpy.zeros((x_size, y_size)) + 10.0
    image = dg.drawGaussians((x_size, y_size),
                             numpy.array([[20.0, 40.0, height, sigma, sigma]]))
    image += background

    # Configure fitter.
    mfit = daoFitC.MultiFitter2D(sigma_range = [1.0, 2.0])
    mfit.initializeC(image)
    mfit.newImage(image)
    mfit.newBackground(background)

    # Configure finder.
    z_values = [0.0]
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2 * sigma,
                                threshold = background[0,0] + 0.5*height,
                                z_values = z_values)
    
    # Find peaks.
    [x, y, z] = mxf.findMaxima([image])

    sigma = numpy.ones(x.size) * sigma
    peaks = {"x" : x,
             "y" : y,
             "z" : z,
             "sigma" : sigma}

    # Pass peaks to fitter.
    mfit.newPeaks(peaks, "finder")

    # Check height.
    h = mfit.getPeakProperty("height")
    for i in range(h.size):
        assert (abs(h[i] - height)/height < 0.1)

    # Check background.
    bg = mfit.getPeakProperty("background")
    for i in range(bg.size):
        assert (abs(bg[i] - 10.0) < 1.0e-6)

    mfit.cleanup(verbose = False)

def test_ia_util_4():
    """
    Multiple z planes test.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64),
              numpy.zeros((x_size,y_size), dtype = numpy.float64),
              numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [1.0,2.0,3.0]

    images[0][20,10] = 1.3
    images[1][20,10] = 1.2
    images[2][20,10] = 1.5

    # Default z range (the entire stack).
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert(x.size == 1)
    assert(abs(z[0] - z_values[2]) < 1.0e-6)

    # z range is limited to adjacent slices.
    #
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_range = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert(x.size == 2)
    assert(abs(z[0] - z_values[0]) < 1.0e-6)
    assert(abs(z[1] - z_values[2]) < 1.0e-6)

    # z range is limited to current slice.
    #
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_range = 0,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert(x.size == 3)
    for i in range(z.size):
        assert(abs(z[i] - z_values[i]) < 1.0e-6)
    
def test_ia_util_5():
    """
    Test that limits on the number of peak duplicates are working.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    # Single peak.
    images[0][10,21] = 1.1
    images[0][10,20] = 1.2
    
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_values = z_values)

    # Find the peak.
    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert (x.size == 1)

    # This should not find anything, since the peak was
    # already found above.
    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert (x.size == 0)

    # Reset and now we should find it again.
    mxf.resetTaken()
    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert (x.size == 1)
    
def test_ia_util_6():
    """
    Test radius.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    images[0][10,23] = 1.2
    images[0][10,22] = 1.1
    images[0][10,21] = 1.1
    images[0][10,20] = 1.2

    # Should only find 1 peak.    
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 3,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert (x.size == 1)

    # Should find two peaks.
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert (x.size == 2)

def test_ia_util_7():
    """
    Test margin.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    images[0][10,5] = 1.2
    images[0][10,6] = 1.1
    images[0][10,7] = 1.1
    images[0][20,10] = 1.1
    images[0][20,11] = 1.2

    # Should only find 1 peak.
    mxf = iaUtilsC.MaximaFinder(margin = 6,
                                radius = 3,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert (x.size == 1)

    # Should find two peaks.
    mxf = iaUtilsC.MaximaFinder(margin = 4,
                                radius = 3,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(imagesCopy(images))
    assert (x.size == 2)

def test_ia_util_8():
    """
    Test allowing multiple duplicates.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    # Single peak.
    images[0][10,21] = 1.1
    images[0][10,20] = 1.2
    
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                n_duplicates = 2,
                                radius = 2,
                                threshold = 1,
                                z_values = z_values)

    # Find the peak.
    np = [1, 1, 0]
    for elt in np:
        [x, y, z] = mxf.findMaxima(imagesCopy(images))
        assert (x.size == elt)

    # Reset.
    mxf.resetTaken()
    
    # Test again.
    np = [1, 1, 0]
    for elt in np:
        [x, y, z] = mxf.findMaxima(imagesCopy(images))
        assert (x.size == elt)

def test_ia_util_9():
    """
    Test runningIfHasNeighbors() function.
    """
    # Test 4 of 5 with new neighbors, one in error state.
    c_x = numpy.array([1.0, 2.0, 3.0, 4.0, 5.0])
    c_y = numpy.array([1.0, 1.0, 1.0, 1.0, 1.0])
    n_x = numpy.array([1.1, 2.1, 3.1, 4.1])
    n_y = numpy.array([1.1, 1.1, 1.1, 1.1])
    status = numpy.array([0, 1, 2, 1, 1], dtype = numpy.int32)

    new_status = iaUtilsC.runningIfHasNeighbors(status, c_x, c_y, n_x, n_y, 0.5)
    correct = [0, 0, 2, 0, 1]
    for i in range(new_status.size):
        assert(new_status[i] == correct[i])

    # Test 2 of 5 with new neighbors, one in error state.
    n_x = numpy.array([1.9, 2.1])
    n_y = numpy.array([1.1, 1.1])
    status = numpy.array([0, 1, 2, 1, 1], dtype = numpy.int32)

    new_status = iaUtilsC.runningIfHasNeighbors(status, c_x, c_y, n_x, n_y, 0.5)
    correct = [0, 0, 2, 1, 1]
    for i in range(new_status.size):
        assert(new_status[i] == correct[i])

    # Test 1 of 2 with new neighbors, but both with radius of each other.
    c_x = numpy.array([2.0, 3.0])
    c_y = numpy.array([2.0, 2.0])
    n_x = numpy.array([1.0])
    n_y = numpy.array([2.0])
    status = numpy.array([1, 1], dtype = numpy.int32)

    new_status = iaUtilsC.runningIfHasNeighbors(status, c_x, c_y, n_x, n_y, 1.5)
    correct = [0, 1]
    for i in range(new_status.size):
        assert(new_status[i] == correct[i])
        
def test_ia_util_10():
    """
    Test markDimmerPeaks() function.
    """
    n_peaks = 25
    x = numpy.random.uniform(size = n_peaks)
    y = numpy.random.uniform(size = n_peaks)
    h = numpy.random.uniform(size = n_peaks) + 1.0
    status = numpy.ones(n_peaks, dtype = numpy.int32)*iaUtilsC.CONVERGED

    # Make first peak the tallest.
    h[0] = 4.0

    # Move last peak outside of the removal radius.
    x[-1] = 4.0

    # Move second to last peak away from everything.
    x[-2] = 40.0
    
    assert (iaUtilsC.markDimmerPeaks(x, y, h, status, 2.0, 5.0) == (n_peaks - 3))

    for i in range(1,n_peaks-2):
        assert(status[i] == iaUtilsC.ERROR)
    assert(status[0] == iaUtilsC.RUNNING)
    assert(status[-1] == iaUtilsC.RUNNING)
    assert(status[-2] == iaUtilsC.CONVERGED)

def test_ia_util_11():
    """
    Test finding peaks in an image.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    # A row of peaks greater than radius with decreasing heights, there
    # should still only be a single maxima.
    images[0][10,11] = 1.6
    images[0][10,12] = 1.5
    images[0][10,13] = 1.4
    images[0][10,14] = 1.3
    images[0][10,15] = 1.2
    images[0][10,16] = 1.1
    
    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z] = mxf.findMaxima(images)
    assert(x.size == 1)
    assert(abs(x[0] - 11.0) < 1.0e-6)
    assert(abs(y[0] - 10.0) < 1.0e-6)

def test_ia_util_12():
    """
    Test removeNeighbors().
    """
    x = numpy.array([1.0, 2.0, 4.0])
    y = numpy.ones(x.size)

    [px, py] = iaUtilsC.removeNeighbors(x, y, 0.5)
    assert(px.size == 3)
    assert(py.size == 3)
    
    [px, py] = iaUtilsC.removeNeighbors(x, y, 1.5)
    assert(px.size == 1)
    assert(py.size == 1)

def test_ia_util_13():
    """
    Test marking peak status (and neighbors status) based on significance.
    """
    px = numpy.array([30.0,35.0,40.0])
    py = numpy.array([30.0,30.0,30.0]) 
    sig = numpy.array([1.0,2.0,3.0])

    # This should not update anything.
    status = numpy.ones(3, dtype = numpy.int32)*iaUtilsC.CONVERGED
    n_marked = iaUtilsC.markLowSignificancePeaks(px,
                                                 py,
                                                 sig,
                                                 status,
                                                 0.0,
                                                 7.0)

    assert(n_marked == 0)
    for i in range(3):
        assert(status[i] == iaUtilsC.CONVERGED)

    # This should only mark the first peak for removal.
    status = numpy.ones(3, dtype = numpy.int32)*iaUtilsC.CONVERGED
    n_marked = iaUtilsC.markLowSignificancePeaks(px,
                                                 py,
                                                 sig,
                                                 status,
                                                 sig[0] + 0.1,
                                                 3.0)
    assert(n_marked == 1)
    assert(status[0] == iaUtilsC.ERROR)
    assert(status[1] == iaUtilsC.CONVERGED)
    assert(status[2] == iaUtilsC.CONVERGED)

    # This should the first peak for removal and the second as RUNNING.
    status = numpy.ones(3, dtype = numpy.int32)*iaUtilsC.CONVERGED
    n_marked = iaUtilsC.markLowSignificancePeaks(px,
                                                 py,
                                                 sig,
                                                 status,
                                                 sig[0] + 0.1,
                                                 7.0)
    assert(n_marked == 1)
    assert(status[0] == iaUtilsC.ERROR)
    assert(status[1] == iaUtilsC.RUNNING)
    assert(status[2] == iaUtilsC.CONVERGED)

    # This should the first peak for removal and both others as RUNNING.
    status = numpy.ones(3, dtype = numpy.int32)*iaUtilsC.CONVERGED
    n_marked = iaUtilsC.markLowSignificancePeaks(px,
                                                 py,
                                                 sig,
                                                 status,
                                                 sig[0] + 0.1,
                                                 11.0)
    assert(n_marked == 1)
    assert(status[0] == iaUtilsC.ERROR)
    assert(status[1] == iaUtilsC.RUNNING)
    assert(status[2] == iaUtilsC.RUNNING)

    # This should the first peak & second for removal and the last as RUNNING.
    status = numpy.ones(3, dtype = numpy.int32)*iaUtilsC.CONVERGED
    n_marked = iaUtilsC.markLowSignificancePeaks(px,
                                                 py,
                                                 sig,
                                                 status,
                                                 sig[1] + 0.1,
                                                 7.0)
    assert(n_marked == 2)
    assert(status[0] == iaUtilsC.ERROR)
    assert(status[1] == iaUtilsC.ERROR)
    assert(status[2] == iaUtilsC.RUNNING)    

def test_ia_util_14():
    """
    Some tests of our KDTree.
    """
    x1 = numpy.array([1.0, 2.0, 3.0])
    y1 = numpy.array([1.0, 1.0, 1.0])
    kd = iaUtilsC.KDTree(x = x1, y = y1)

    x2 = numpy.array([1.1, 4.0])
    y2 = numpy.array([1.0, 1.0])

    [dist, index] = kd.nearest(x2, y2, 0.2)
    assert(abs(dist[0] - 0.1) < 1.0e-6)
    assert(abs(dist[1] + 1.0) < 1.0e-6)
    assert(index[0] == 0)
    assert(index[1] == -1)
    
    [dist, index] = kd.nearest(x2, y2, 3.1)
    assert(abs(dist[0] - 0.1) < 1.0e-6)
    assert(abs(dist[1] - 1.0) < 1.0e-6)
    assert(index[0] == 0)
    assert(index[1] == 2)    

    kd.cleanup()

def test_ia_util_15():
    """
    Test non-integer radius values.
    """
    x_size = 100
    y_size = 80
    images = [numpy.zeros((x_size,y_size), dtype = numpy.float64)]
    z_values = [0.1]

    images[0][8,10] = 1.1
    images[0][10,10] = 1.5
    images[0][11,11] = 1.1

    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 1.0,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z, h] = mxf.findMaxima(images, want_height = True)
    assert(numpy.allclose(x, numpy.array([10,10,11])))
    assert(numpy.allclose(y, numpy.array([8,10,11])))

    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 1.99,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z, h] = mxf.findMaxima(images, want_height = True)
    assert(numpy.allclose(x, numpy.array([10,10])))
    assert(numpy.allclose(y, numpy.array([8,10])))

    mxf = iaUtilsC.MaximaFinder(margin = 1,
                                radius = 2.01,
                                threshold = 1,
                                z_values = z_values)

    [x, y, z, h] = mxf.findMaxima(images, want_height = True)
    assert(numpy.allclose(x, numpy.array([10])))
    assert(numpy.allclose(y, numpy.array([10])))


if (__name__ == "__main__"):
    test_ia_util_1()
    test_ia_util_2()
    test_ia_util_3()
    test_ia_util_4()
    test_ia_util_5()
    test_ia_util_6()
    test_ia_util_7()
    test_ia_util_8()
    test_ia_util_9()
    test_ia_util_10()
    test_ia_util_11()
    test_ia_util_12()
    test_ia_util_13()
    test_ia_util_14()
    test_ia_util_15()
