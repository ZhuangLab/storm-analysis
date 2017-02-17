#!/usr/bin/env python

import storm_analysis

import storm_analysis.test.verifications as veri

    
def test_setup_A_matrix():

    # Test setupAMatrix.
    a_matrix_file = storm_analysis.getPathOutputTest("test_l1h")
    storm_analysis.removeFile(a_matrix_file)

    from storm_analysis.L1H.setup_A_matrix import setupAMatrix
    setupAMatrix("theoritical", a_matrix_file, 1.0, False)


def test_l1h():
    
    # Test L1H.
    movie_name = storm_analysis.getData("test/data/test_l1h.dax")
    settings = storm_analysis.getData("test/data/test_l1h.xml")
    hres = storm_analysis.getPathOutputTest("test_l1h_list.hres")
    mlist = storm_analysis.getPathOutputTest("test_l1h_list.bin")

    storm_analysis.removeFile(hres)
    storm_analysis.removeFile(mlist)

    from storm_analysis.L1H.cs_analysis import analyze
    analyze(movie_name, settings, hres, mlist)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if (num_locs != 1986):
        raise Exception("L1H did not find the expected number of localizations.")
    

if (__name__ == "__main__"):
    test_setup_A_matrix()
    test_l1h()
