#!/usr/bin/env python

import storm_analysis

def test_measure_psf():

    movie = storm_analysis.getData("test/data/test_spliner.dax")
    mlist = storm_analysis.getData("test/data/test_spliner_olist.bin")
    psf = storm_analysis.getPathOutputTest("test_spliner_psf.psf")
    storm_analysis.removeFile(psf)

    from storm_analysis.spliner.measure_psf import measurePSF
    measurePSF(movie, "", mlist, psf)

    
def test_psf_to_spline():

    psf = storm_analysis.getPathOutputTest("test_spliner_psf.psf")
    spline = storm_analysis.getPathOutputTest("test_spliner_psf.spline")
    storm_analysis.removeFile(spline)
    
    from storm_analysis.spliner.psf_to_spline import psfToSpline
    psfToSpline(psf, spline, 10)


def test_spliner_std():

    movie_name = storm_analysis.getData("test/data/test_spliner.dax")
    settings = storm_analysis.getData("test/data/test_spliner_dh.xml")
    mlist = storm_analysis.getPathOutputTest("test_spliner_slist.bin")
    storm_analysis.removeFile(mlist)
        
    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)

    
def test_spliner_std():

    movie_name = storm_analysis.getData("test/data/test_300x200_dh.dax")
    settings = storm_analysis.getData("test/data/test_spliner_dh.xml")
    mlist = storm_analysis.getPathOutputTest("test_300x200_dh_slist.bin")
    storm_analysis.removeFile(mlist)

    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)
    

def test_spliner_fista():

    movie_name = storm_analysis.getData("test/data/test_spliner.dax")
    settings = storm_analysis.getData("test/data/test_spliner_dh_fista.xml")
    mlist = storm_analysis.getPathOutputTest("test_spliner_flist.bin")
    storm_analysis.removeFile(mlist)

    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)


def test_spliner_fista_non_square():

    movie_name = storm_analysis.getData("test/data/test_300x200_dh.dax")
    settings = storm_analysis.getData("test/data/test_spliner_dh_fista.xml")
    mlist = storm_analysis.getPathOutputTest("test_300x200_dh_flist.bin")
    storm_analysis.removeFile(mlist)

    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)

    
if (__name__ == "__main__"):
    test_measure_psf()
    test_psf_to_spline()
    test_spliner_std()
    test_spliner_std_non_square()
    test_spliner_fista()
    test_spliner_fista_non_square()

