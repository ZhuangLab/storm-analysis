#!/usr/bin/env python

import storm_analysis


def test_wavelet_bgr():

    movie_in = storm_analysis.getData("test/data/test_bg_sub.dax")
    movie_out = storm_analysis.getPathOutputTest("test_bg_sub_wbgr.dax")

    from storm_analysis.wavelet_bgr.wavelet_bgr import waveletBGRSub
    waveletBGRSub(movie_in, movie_out, "db4", 2, 2, 10)


if (__name__ == "__main__"):
    test_wavelet_bgr()


