#!/usr/bin/env python
"""
Tests for the edge_fill parameter, how the padding margin around each
frame is filled for Spliner.
"""
import numpy
import pytest
import scipy.ndimage

import storm_analysis

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.spliner.find_peaks_std as findPeaksSTD

import storm_analysis.test.verifications as veri


def splinerParameters(edge_fill = None, edge_fill_weight = None):
    p = params.ParametersSplinerSTD().initFromFile(storm_analysis.getData("test/data/test_spliner_dh.xml"))
    if edge_fill is not None:
        p.setAttr("edge_fill", "string", edge_fill)
    if edge_fill_weight is not None:
        p.setAttr("edge_fill_weight", "float", edge_fill_weight)
    return p


def test_edge_pad_array():
    image = numpy.arange(12, dtype = numpy.float64).reshape(3, 4)
    padded = fitting.edgePadArray(image, 2)

    assert (padded.shape == (7, 8))
    assert numpy.array_equal(padded[2:5,2:6], image)

    # Each side repeats its edge pixels, the corners the corner pixels.
    for i in range(2):
        assert numpy.array_equal(padded[i,2:6], image[0,:])
        assert numpy.array_equal(padded[6-i,2:6], image[-1,:])
        assert numpy.array_equal(padded[2:5,i], image[:,0])
        assert numpy.array_equal(padded[2:5,7-i], image[:,-1])
    assert numpy.all(padded[:2,:2] == image[0,0])
    assert numpy.all(padded[5:,6:] == image[-1,-1])

    # No margin, no change.
    assert (fitting.edgePadArray(image, 0) is image)


def test_planar_background_plane():
    """
    The planar estimate extends a plane exactly, in the frame and in the margin.
    """
    for margin in [10, 18, 32]:
        shape = (64 + 2*margin, 80 + 2*margin)
        [rows, columns] = numpy.indices(shape)
        plane = 100.0 + 0.5*rows - 0.25*columns

        # Only the real pixels count, so what is in the margin doesn't matter.
        image = numpy.zeros(shape)
        image[margin:-margin,margin:-margin] = plane[margin:-margin,margin:-margin]

        planar = fitting.PlanarBackground(shape = shape, margin = margin, sigma = 8.0)
        estimate = planar.estimate(image)
        valid = planar.valid
        planar.cleanup()

        # Exact everywhere a fit window can reach, 16 pixels into the margin
        # for Spliner. Farther out, where only a few rows of real pixels are
        # within the Gaussian's reach, the fit amplifies round-off (to about
        # 1e-4 at 31 pixels into a 32 pixel margin), but nothing is fitted
        # there.
        reach = max(margin - 16, 0)
        reachable = (slice(reach, shape[0] - reach), slice(reach, shape[1] - reach))
        assert numpy.all(valid[reachable])
        assert numpy.allclose(estimate[reachable], plane[reachable], rtol = 0.0, atol = 1.0e-8)

        # With a 32 pixel margin the outermost row is at the Gaussian's 4
        # sigma reach, a plane can't be fitted there, and it gets the median.
        if (margin == 32):
            assert not numpy.any(valid[0,:])
            assert numpy.allclose(estimate[0,:], numpy.median(image[margin:-margin,margin:-margin]))


def test_planar_background_interior():
    """
    Where the Gaussian lies entirely on real pixels the planar estimate is
    a Gaussian low pass filter.
    """
    margin = 18
    shape = (100 + 2*margin, 120 + 2*margin)
    image = numpy.random.default_rng(1).normal(100.0, 10.0, shape)

    planar = fitting.PlanarBackground(shape = shape, margin = margin, sigma = 8.0)
    estimate = planar.estimate(image)
    planar.cleanup()

    lowpass = scipy.ndimage.gaussian_filter(image, 8.0, truncate = 4.0)
    interior = (slice(margin + 32, -(margin + 32)), slice(margin + 32, -(margin + 32)))
    assert numpy.allclose(estimate[interior], lowpass[interior], rtol = 0.0, atol = 1.0e-9)


def test_edge_fill_parameter_checks():
    with pytest.raises(fitting.FittingException):
        findPeaksSTD.initFindAndFit(splinerParameters(edge_fill = "constant"))

    with pytest.raises(fitting.FittingException):
        findPeaksSTD.initFindAndFit(splinerParameters(edge_fill = "edge", edge_fill_weight = 0.5))

    # The default weight is 30.
    finder = findPeaksSTD.initFindAndFit(splinerParameters(edge_fill = "planar"))
    assert (finder.peak_finder.edge_fill_weight == 30.0)
    finder.cleanUp()

    # And the default fill is the mirror.
    finder = findPeaksSTD.initFindAndFit(splinerParameters())
    assert (finder.peak_finder.edge_fill == "mirror")
    finder.cleanUp()


def test_edge_fill_margin_variance():
    """
    With the 'planar' and 'edge' fills the C fitter's variance term in the
    margin is (edge_fill_weight - 1) times the fill, frame by frame.
    """
    for edge_fill in ["planar", "edge"]:
        p = splinerParameters(edge_fill = edge_fill, edge_fill_weight = 10.0)
        finder = findPeaksSTD.initFindAndFit(p)
        frame_reader = analysisIO.FrameReaderStd(movie_file = storm_analysis.getData("test/data/test_spliner.dax"),
                                                 parameters = p)
        movie_reader = analysisIO.MovieReader(frame_reader = frame_reader, parameters = p)
        movie_reader.setup(0)
        movie_reader.max_frame = 3

        while movie_reader.nextFrame():
            finder.analyzeImage(movie_reader)

            image = finder.peak_finder.image
            margin = finder.peak_finder.margin
            mfitter = finder.peak_fitter.mfitter
            scmos_term = numpy.ctypeslib.as_array(mfitter.mfit.contents.scmos_term,
                                                  shape = (image.size,)).reshape(image.shape)

            expected = numpy.maximum(9.0*image, 0.0)
            expected[margin:-margin,margin:-margin] = 0.0
            assert numpy.allclose(scmos_term, expected)
            assert numpy.any(scmos_term > 0.0)

        finder.cleanUp()


def analyzeWithEdgeFill(edge_fill):
    """
    Run Spliner on the test movie with an edge fill, return the localizations.
    """
    name = "test_edge_fill_" + str(edge_fill)
    xml = storm_analysis.getPathOutputTest(name + ".xml")
    mlist = storm_analysis.getPathOutputTest(name + ".hdf5")
    storm_analysis.removeFile(mlist)

    p = splinerParameters(edge_fill = edge_fill)
    p.toXMLFile(xml, remove_paths = False)

    from storm_analysis.spliner.spline_analysis import analyze
    analyze(storm_analysis.getData("test/data/test_spliner.dax"), mlist, xml)

    locs = []
    with saH5Py.SAH5Py(mlist) as h5:
        assert h5.isAnalysisFinished()
        for fnum, frame_locs in h5.localizationsIterator(fields = ["x", "y", "z"]):
            locs.append(frame_locs)
    return locs


def test_edge_fill_spliner():
    """
    The mirror is the default, and the other fills find about as many
    localizations on a movie where few emitters are near the edge.
    """
    default = analyzeWithEdgeFill(None)
    mirror = analyzeWithEdgeFill("mirror")
    assert (len(default) == len(mirror))
    for [locs1, locs2] in zip(default, mirror):
        for field in locs1:
            assert numpy.array_equal(locs1[field], locs2[field])

    n_mirror = sum([frame_locs["x"].size for frame_locs in mirror])
    for edge_fill in ["planar", "edge"]:
        locs = analyzeWithEdgeFill(edge_fill)
        n_locs = sum([frame_locs["x"].size for frame_locs in locs])
        assert veri.verifyIsCloseEnough(n_locs, n_mirror, margin = 0.02)


if (__name__ == "__main__"):
    test_edge_pad_array()
    test_planar_background_plane()
    test_planar_background_interior()
    test_edge_fill_parameter_checks()
    test_edge_fill_margin_variance()
    test_edge_fill_spliner()
