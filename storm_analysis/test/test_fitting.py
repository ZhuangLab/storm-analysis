#!/usr/bin/env python
import numpy
import os

import storm_analysis

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.parameters as params


def test_fitting_exception():
    """
    Test that FittingException is actually an exception.

    It was declared with 'def' rather than 'class', which made it a function
    and not an exception at all.
    """
    assert (isinstance(fitting.FittingException, type))
    assert (issubclass(fitting.FittingException, Exception))


def test_fitting_no_peak_finder():
    """
    A PeakFinder sub-class that does not provide peakFinder() should get
    FittingException, and not a TypeError about the exception not deriving
    from BaseException. The method does not use 'self', so this can be
    checked without building a complete finder.
    """
    try:
        fitting.PeakFinder.peakFinder(None, None)
    except fitting.FittingException:
        pass
    else:
        assert(False)


def test_fitting_peak_mask_circular():
    """
    A circular AOI on a non-square image.

    Both coordinate ranges were built from shape[0], so on a non-square image
    the mask did not match the array it was indexing.
    """
    margin = 10
    [xc, yc, radius] = [150, 40, 20]

    p = params.ParametersDAO()
    p.setAttr("x_center", "int", xc)
    p.setAttr("y_center", "int", yc)
    p.setAttr("aoi_radius", "int", radius)

    # More columns than rows, so that the two axes cannot be confused.
    shape = (100, 200)
    mask = fitting.peakMask(shape, p, margin)

    assert (mask.shape == shape)

    [rows, cols] = numpy.nonzero(mask)

    # x is the column axis and y is the row axis, as in the square AOI branch.
    assert (abs(numpy.mean(rows) - (yc + margin)) < 0.1)
    assert (abs(numpy.mean(cols) - (xc + margin)) < 0.1)

    # And it should be a disc of the requested radius.
    assert (abs(rows.size - numpy.pi*radius*radius) < 0.02*numpy.pi*radius*radius)


def test_fitting_peak_mask_square_aoi():
    """
    The rectangular AOI branch masks x by column and y by row.
    """
    margin = 5
    p = params.ParametersDAO()
    p.setAttr("x_start", "int", 20)
    p.setAttr("y_start", "int", 10)

    mask = fitting.peakMask((60, 120), p, margin)

    # Everything left of x_start + margin is masked, and nothing beyond it is.
    assert (numpy.count_nonzero(mask[:,:20+margin]) == 0)
    assert (numpy.count_nonzero(mask[10+margin:,20+margin:]) > 0)

    # Everything above y_start + margin is masked.
    assert (numpy.count_nonzero(mask[:10+margin,:]) == 0)


def test_fitting_check_mode():
    """
    check_mode writes diagnostic images from the arbitrary PSF peak finder,
    which is the one Spliner, Pupilfn and PSF FFT all use. It called
    tifffile.imsave() and TiffWriter.save(), both of which have been removed
    from tifffile, so turning it on raised AttributeError.
    """
    import storm_analysis.spliner.find_peaks_std as findPeaksStd

    settings = storm_analysis.getData("test/data/test_spliner_2D.xml")
    parameters = params.ParametersSplinerSTD().initFromFile(settings, warnings = False)
    parameters.setAttr("spline", "filename",
                       storm_analysis.getData("test/data/test_spliner_psf_2d.spline"))

    find_fit = findPeaksStd.initFindAndFit(parameters)
    finder = find_fit.peak_finder

    image = numpy.random.normal(100.0, 5.0, (100, 100))
    image = finder.padArray(image)

    # check_mode writes into the working directory.
    cwd = os.getcwd()
    os.chdir(storm_analysis.getPathOutputTest(""))
    try:
        finder.check_mode = True
        finder.newImage(image)
        fit_peaks_image = numpy.zeros(image.shape)
        finder.estimateBackground(fit_peaks_image, None)

        # This is what used to raise AttributeError.
        finder.peakFinder(fit_peaks_image)
    finally:
        os.chdir(cwd)
        find_fit.cleanUp()


if (__name__ == "__main__"):
    test_fitting_exception()
    test_fitting_no_peak_finder()
    test_fitting_peak_mask_circular()
    test_fitting_peak_mask_square_aoi()
    test_fitting_check_mode()
