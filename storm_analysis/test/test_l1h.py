#!/usr/bin/env python

import storm_analysis


def test_setup_A_matrix():

    a_matrix_file = storm_analysis.get_data("test/data/test_l1h")

    from storm_analysis.L1H.setup_A_matrix import setupAMatrix
    setupAMatrix("theoritical", a_matrix_file, 1.0, False)

    
def test_l1h():

    movie_name = storm_analysis.get_data("test/data/test_l1h.dax")
    settings = storm_analysis.get_data("test/data/test_l1h.xml")
    hres = storm_analysis.get_path_output_test("test_l1h_list.hres")
    mlist = storm_analysis.get_path_output_test("test_l1h_list.bin")

    print(settings)
    
    from storm_analysis.L1H.cs_analysis import analyze
    analyze(movie_name, settings, hres, mlist)
    

if (__name__ == "__main__"):
    test_setup_A_matrix()
    test_l1h()

