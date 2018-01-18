#!/usr/bin/env python
"""
Tests for Spliner analysis.
"""
import sys

import storm_analysis

import storm_analysis.test.verifications as veri


def test_measure_psf():

    movie = storm_analysis.getData("test/data/test_spliner.dax")
    mlist = storm_analysis.getData("test/data/test_spliner_ref.hdf5")
    psf = storm_analysis.getPathOutputTest("test_spliner_psf.psf")
    storm_analysis.removeFile(psf)

    from storm_analysis.spliner.measure_psf import measurePSF
    measurePSF(movie, "", mlist, psf)

    
def test_measure_psf_2D():

    movie = storm_analysis.getData("test/data/test.dax")
    mlist = storm_analysis.getData("test/data/test_ref.hdf5")
    psf = storm_analysis.getPathOutputTest("test_spliner_psf_2d.psf")
    storm_analysis.removeFile(psf)

    from storm_analysis.spliner.measure_psf import measurePSF
    measurePSF(movie, "", mlist, psf, want2d = True, aoi_size = 5)

    
def _test_psf_to_spline():

    psf = storm_analysis.getPathOutputTest("test_spliner_psf.psf")
    spline = storm_analysis.getPathOutputTest("test_spliner_psf.spline")
    storm_analysis.removeFile(spline)
    
    from storm_analysis.spliner.psf_to_spline import psfToSpline
    psfToSpline(psf, spline, 10)

    
def _test_psf_to_spline_2D():

    psf = storm_analysis.getPathOutputTest("test_spliner_psf_2d.psf")
    spline = storm_analysis.getPathOutputTest("test_spliner_psf_2d.spline")
    storm_analysis.removeFile(spline)
    
    from storm_analysis.spliner.psf_to_spline import psfToSpline
    psfToSpline(psf, spline, 7)

    
def test_spliner_std():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
                    
    movie_name = storm_analysis.getData("test/data/test_spliner.dax")
    settings = storm_analysis.getData("test/data/test_spliner_dh.xml")
    mlist = storm_analysis.getPathOutputTest("test_spliner_dh.hdf5")
    storm_analysis.removeFile(mlist)
        
    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 720):
        raise Exception("Spliner 3D did not find the expected number of localizations.")        


def test_spliner_std_2D():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_spliner_2D.xml")
    mlist = storm_analysis.getPathOutputTest("test_spliner_2D.hdf5")
    storm_analysis.removeFile(mlist)
        
    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)
    
    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 2004):
        raise Exception("Spliner 2D did not find the expected number of localizations.")

    
def test_spliner_std_non_square():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    movie_name = storm_analysis.getData("test/data/test_300x200_dh.dax")
    settings = storm_analysis.getData("test/data/test_spliner_dh.xml")
    mlist = storm_analysis.getPathOutputTest("test_spliner_dh.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 120):               
        raise Exception("Spliner 3D non square did not find the expected number of localizations.")
    

def _test_spliner_fista():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return    

    movie_name = storm_analysis.getData("test/data/test_spliner.dax")
    settings = storm_analysis.getData("test/data/test_spliner_dh_fista.xml")
    mlist = storm_analysis.getPathOutputTest("test_spliner_dh_fista.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 36):               
        raise Exception("Spliner 3D FISTA did not find the expected number of localizations.")


def _test_spliner_fista_2D():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_spliner_2D_fista.xml")
    mlist = storm_analysis.getPathOutputTest("test_spliner_2D_fista.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 587):
        raise Exception("Spliner 2D FISTA did not find the expected number of localizations.")


def _test_spliner_fista_non_square():

    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    movie_name = storm_analysis.getData("test/data/test_300x200_dh.dax")
    settings = storm_analysis.getData("test/data/test_spliner_dh_fista.xml")
    mlist = storm_analysis.getPathOutputTest("test_spliner_dh_fista.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.spliner.spline_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 24):               
        raise Exception("Spliner 3D FISTA non square did not find the expected number of localizations.")
    
    
if (__name__ == "__main__"):
    test_measure_psf()
    test_measure_psf_2D()
#    _test_psf_to_spline()
#    _test_psf_to_spline_2D()
    test_spliner_std()
    test_spliner_std_2D()
    test_spliner_std_non_square()
    _test_spliner_fista()
    _test_spliner_fista_2D()
    _test_spliner_fista_non_square()
