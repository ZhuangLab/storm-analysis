#!/usr/bin/env python

import storm_analysis

import storm_analysis.sa_library.readinsight3 as readinsight3

import storm_analysis.test.verifications as veri

    
def test_homotopy_psf():

    movie = storm_analysis.getData("test/data/test.dax")
    mlist = storm_analysis.getData("test/data/test_olist.bin")
    psf = storm_analysis.getPathOutputTest("l1h_psf.psf")
    storm_analysis.removeFile(psf)

    from storm_analysis.L1H.homotopy_psf import homotopyPSF
    homotopyPSF(movie, mlist, psf)
    

def test_setup_A_matrix():

    # Test setupAMatrix.
    a_matrix_file = storm_analysis.getPathOutputTest("test_l1h")
    storm_analysis.removeFile(a_matrix_file)

    from storm_analysis.L1H.setup_A_matrix import setupAMatrix
    setupAMatrix("theoritical", a_matrix_file, 1.0, False)


def _test_l1h():
    
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
    #
    # FIXME: Change L1H to use the HDF5 format.
    #
    num_locs = readinsight3.loadI3File(mlist)["x"].size
    if not veri.verifyIsCloseEnough(num_locs, 1986):        
        raise Exception("L1H did not find the expected number of localizations.")
    

if (__name__ == "__main__"):
    test_homotopy_psf()
    test_setup_A_matrix()
    _test_l1h()
    
