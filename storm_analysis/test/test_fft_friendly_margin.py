#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.sa_h5py as saH5Py


def test_max_prime_factor():
    assert (fitting.maxPrimeFactor(1) == 1)
    assert (fitting.maxPrimeFactor(2) == 2)
    assert (fitting.maxPrimeFactor(256) == 2)
    assert (fitting.maxPrimeFactor(270) == 5)
    assert (fitting.maxPrimeFactor(280) == 7)

    # The size a 256x256 movie pads to by default, and the reason this
    # module exists. 268 = 2*2*67.
    assert (fitting.maxPrimeFactor(268) == 67)


def test_fft_friendly_margin():

    # Common sensor sizes should all move off a large prime factor.
    for size in [256, 512, 1024, 2048]:
        margin = fitting.fftFriendlyMargin((size, size), 6)
        assert (margin >= 6)
        assert (fitting.maxPrimeFactor(size + 2*margin) <= 7)

    # 268 is the case that motivated this, 270 is the answer.
    assert (fitting.fftFriendlyMargin((256, 256), 6) == 7)

    # Non-square images have to satisfy both axes with one margin.
    margin = fitting.fftFriendlyMargin((256, 300), 6)
    assert (fitting.maxPrimeFactor(256 + 2*margin) <= 7)
    assert (fitting.maxPrimeFactor(300 + 2*margin) <= 7)

    # A margin that is already fine is left alone.
    assert (fitting.fftFriendlyMargin((256, 256), 7) == 7)

    # If nothing suitable is in range the original margin comes back, so
    # this can decline to help but never fail.
    assert (fitting.fftFriendlyMargin((256, 256), 6, max_prime = 2,
                                      max_extra = 2) == 6)


def test_fft_friendly_margin_invariant():
    """
    Growing the margin must not change the localizations, only how quickly
    they are found.
    """
    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")

    from storm_analysis.daostorm_3d.mufit_analysis import analyze

    results = []
    for enabled in [1, 0]:
        xml = storm_analysis.getPathOutputTest("test_ffm_" + str(enabled) + ".xml")
        mlist = storm_analysis.getPathOutputTest("test_ffm_" + str(enabled) + ".hdf5")
        storm_analysis.removeFile(mlist)

        p = params.ParametersDAO().initFromFile(settings)
        p.setAttr("fft_friendly_margin", "int", enabled)
        p.toXMLFile(xml)

        analyze(movie_name, mlist, xml)

        with saH5Py.SAH5Py(mlist) as h5:
            locs = h5.getLocalizations()
        results.append(locs)

    [with_ffm, without_ffm] = results

    assert (with_ffm["x"].size == without_ffm["x"].size)

    # Different FFT sizes reorder the floating point work, so these agree to
    # round-off rather than exactly. The tolerance is still many orders of
    # magnitude below any meaningful localization precision.
    for field in ["x", "y"]:
        assert (numpy.allclose(with_ffm[field], without_ffm[field],
                               atol = 1.0e-3, rtol = 0.0)), field


if (__name__ == "__main__"):
    test_max_prime_factor()
    test_fft_friendly_margin()
    test_fft_friendly_margin_invariant()
