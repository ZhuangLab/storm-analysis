#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.daostorm_3d.find_peaks as findPeaks
import storm_analysis.sa_library.parameters as params

import storm_analysis.test.verifications as veri


def test_3ddao_2d_fixed():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_fixed.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1998):
        raise Exception("3D-DAOSTORM 2D fixed did not find the expected number of localizations.")
    

def test_3ddao_2d_fixed_gt():
    """
    Start fitting from ground truth locations.
    """
    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed_gt.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_fixed_gt.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 200):
        raise Exception("3D-DAOSTORM 2D fixed ground truth did not find the expected number of localizations.")


def test_3ddao_2d_fixed_gt_text():
    """
    Start fitting from ground truth locations (text file version).
    """
    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed_gt_text.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_fixed_gt_text.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 200):
        raise Exception("3D-DAOSTORM 2D fixed ground truth did not find the expected number of localizations.")    
    
    
def test_3ddao_2d_fixed_low_snr():

    movie_name = storm_analysis.getData("test/data/test_low_snr.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed_low_snr.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_fixed_low_snr.hdf5")
    storm_analysis.removeFile(mlist)
    
    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 392):
        raise Exception("3D-DAOSTORM 2D fixed low snr did not find the expected number of localizations.")
    

def test_3ddao_2d_fixed_non_square():

    movie_name = storm_analysis.getData("test/data/test_300x200.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_300x200.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1002):
        raise Exception("3D-DAOSTORM 2D fixed non square did not find the expected number of localizations.")
    
    
def test_3ddao_2d():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1970):
        raise Exception("3D-DAOSTORM 2D did not find the expected number of localizations.")


def test_3ddao_3d():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_3d.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_3d.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1956):
        raise Exception("3D-DAOSTORM 3D did not find the expected number of localizations.")    


def test_3ddao_Z():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_Z.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_Z.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1955):
        raise Exception("3D-DAOSTORM Z did not find the expected number of localizations.")


def test_3ddao_scmos_cal():
    """
    Test that scmos calibration data is initialized to 0.0.
    """
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    parameters = params.ParametersDAO().initFromFile(settings)

    # Create analysis object and reach deep into it..
    find_fit = findPeaks.initFindAndFit(parameters)
    fitter = find_fit.peak_fitter
    mfitter = fitter.mfitter

    # Initialize with an image.
    image = numpy.ones((100,100))
    fitter.newImage(image)

    # Verify that the image is still all ones.
    resp = mfitter.getResidual()
    assert(numpy.max(resp - 1.0) < 1.0e-6)

    # Cleanup.
    fitter.cleanUp()

    
if (__name__ == "__main__"):
    test_3ddao_2d_fixed()
    test_3ddao_2d_fixed_gt()
    test_3ddao_2d_fixed_gt_text()
    test_3ddao_2d_fixed_low_snr()
    test_3ddao_2d_fixed_non_square()
    test_3ddao_2d()
    test_3ddao_3d()
    test_3ddao_Z()
    test_3ddao_scmos_cal()

