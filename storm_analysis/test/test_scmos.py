#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.daostorm_3d.find_peaks as findPeaks
import storm_analysis.sa_library.parameters as params

import storm_analysis.test.verifications as veri


def test_scmos_2d_fixed():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_sc_2d_fixed.xml")
    mlist = storm_analysis.getPathOutputTest("test_sc_2d_fixed.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.sCMOS.scmos_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1992):        
        raise Exception("sCMOS 2D fixed did not find the expected number of localizations.")    

    
def test_scmos_2d():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_sc_2d.xml")
    mlist = storm_analysis.getPathOutputTest("test_sc_2d.hdf5")
    storm_analysis.removeFile(mlist)
    
    from storm_analysis.sCMOS.scmos_analysis import analyze    
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1961):
        raise Exception("sCMOS 2D did not find the expected number of localizations.")
    
    
def test_scmos_3d():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_sc_3d.xml")
    mlist = storm_analysis.getPathOutputTest("test_sc_3d.hdf5")
    storm_analysis.removeFile(mlist)
    
    from storm_analysis.sCMOS.scmos_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1950):
        raise Exception("sCMOS 3D did not find the expected number of localizations.")

    # Verify that the Z values actually got calculated.
    if not veri.verifyZWasCalculated(mlist):
        raise Exception("Z values were not calculated for sCMOS 3D fitting.")

def test_scmos_Z():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_sc_Z.xml")
    mlist = storm_analysis.getPathOutputTest("test_sc_Z.hdf5")
    storm_analysis.removeFile(mlist)
    
    from storm_analysis.sCMOS.scmos_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1942):
        raise Exception("sCMOS Z did not find the expected number of localizations.")
    
def test_scmos_scmos_cal():
    """
    Test that scmos calibration data is initialized correctly.
    """
    settings = storm_analysis.getData("test/data/test_sc_2d_fixed.xml")
    parameters = params.ParametersSCMOS().initFromFile(settings)

    # Create analysis object and reach deep into it..
    find_fit = findPeaks.initFindAndFit(parameters)
    fitter = find_fit.peak_fitter
    mfitter = fitter.mfitter

    # Get sCMOS calibration value.
    scmos_cal_value = mfitter.scmos_cal[0,0]
    
    # Initialize with an image.
    image = numpy.ones(mfitter.scmos_cal.shape)
    fitter.newImage(image)

    # Verify that the image has the sCMOS term added.
    resp = mfitter.getResidual()
    assert(numpy.max(resp - (1.0 + scmos_cal_value)) < 1.0e-6)

    # Cleanup.
    fitter.cleanUp()


if (__name__ == "__main__"):
    test_scmos_2d_fixed()
    test_scmos_2d()
    test_scmos_3d()
    test_scmos_Z()
    test_scmos_scmos_cal()
