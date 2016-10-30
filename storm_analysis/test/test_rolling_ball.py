#!/usr/bin/env python

import storm_analysis


def test_rolling_ball():

    movie_in = storm_analysis.getData("test/data/test_bg_sub.dax")
    movie_out = storm_analysis.getPathOutputTest("test_bg_sub_rb.dax")
    
    from storm_analysis.rolling_ball_bgr.rolling_ball import rollingBallSub
    rollingBallSub(movie_in, movie_out, 10, 1)


if (__name__ == "__main__"):
    test_rolling_ball()


