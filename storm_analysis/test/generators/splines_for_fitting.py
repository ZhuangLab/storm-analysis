#!/usr/bin/env python
"""
Creates the splines that Spliner() will use in testing. Once they
are created you'll need to copy them to the test/data directory.

Hazen 09/18
"""
import sys

import storm_analysis
import storm_analysis.spliner.measure_psf as measurePSF
import storm_analysis.spliner.psf_to_spline as psfToSpline
    

def create2DSpline():
    movie = storm_analysis.getData("test/data/test.dax")
    mlist = storm_analysis.getData("test/data/test_ref.hdf5")
    psf = storm_analysis.getPathOutputTest("test_spliner_psf_2d.psf")
    spline = storm_analysis.getPathOutputTest("test_spliner_psf_2d.spline")
    
    storm_analysis.removeFile(psf)
    storm_analysis.removeFile(spline)

    measurePSF.measurePSF(movie, "", mlist, psf, want2d = True, aoi_size = 5)
    psfToSpline.psfToSpline(psf, spline, 4)
    
    
def create3DSpline():

    movie = storm_analysis.getData("test/data/test_spliner.dax")
    mlist = storm_analysis.getData("test/data/test_spliner_ref.hdf5")
    psf = storm_analysis.getPathOutputTest("test_spliner_psf.psf")
    spline = storm_analysis.getPathOutputTest("test_spliner_psf.spline")

    storm_analysis.removeFile(psf)
    storm_analysis.removeFile(spline)
    
    measurePSF.measurePSF(movie, "", mlist, psf, aoi_size = 6)
    psfToSpline.psfToSpline(psf, spline, 5)
    

if (__name__ == "__main__"):
    create2DSpline()
    create3DSpline()
    
