#!/usr/bin/env python
import storm_analysis.sa_library.fitting as fitting


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


if (__name__ == "__main__"):
    test_fitting_exception()
    test_fitting_no_peak_finder()
