#!/usr/bin/env python

import storm_analysis


def test_scmos_2d_fixed():

    movie_name = storm_analysis.get_data("test/data/test.dax")
    settings = storm_analysis.get_data("test/data/test_sc_2d_fixed.xml")
    mlist = storm_analysis.get_path_output_test("test_sc_2d_fixed.bin")

    from storm_analysis.sCMOS.scmos_analysis import analyze
    analyze(movie_name, mlist, settings)

    
def test_scmos_2d():

    movie_name = storm_analysis.get_data("test/data/test.dax")
    settings = storm_analysis.get_data("test/data/test_sc_2d.xml")
    mlist = storm_analysis.get_path_output_test("test_sc_2d.bin")

    from storm_analysis.sCMOS.scmos_analysis import analyze    
    analyze(movie_name, mlist, settings)

    
def test_scmos_3d():

    movie_name = storm_analysis.get_data("test/data/test.dax")
    settings = storm_analysis.get_data("test/data/test_sc_3d.xml")
    mlist = storm_analysis.get_path_output_test("test_sc_3d.bin")

    from storm_analysis.sCMOS.scmos_analysis import analyze
    analyze(movie_name, mlist, settings)


def test_scmos_Z():

    movie_name = storm_analysis.get_data("test/data/test.dax")
    settings = storm_analysis.get_data("test/data/test_sc_Z.xml")
    mlist = storm_analysis.get_path_output_test("test_sc_Z.bin")

    from storm_analysis.sCMOS.scmos_analysis import analyze
    analyze(movie_name, mlist, settings)


if (__name__ == "__main__"):
    test_scmos_2d_fixed()
    test_scmos_2d()
    test_scmos_3d()
    test_scmos_Z()
