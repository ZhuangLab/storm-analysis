#!/usr/bin/env python

import shutil
import storm_analysis
    
def test_l1h():

    # Test setupAMatrix.
    a_matrix_file = storm_analysis.getPathOutputTest("test_l1h")
    storm_analysis.removeFile(a_matrix_file)

    from storm_analysis.L1H.setup_A_matrix import setupAMatrix
    setupAMatrix("theoritical", a_matrix_file, 1.0, False)

    # Test L1H.
    movie_name = storm_analysis.getData("test/data/test_l1h.dax")
    hres = storm_analysis.getPathOutputTest("test_l1h_list.hres")
    mlist = storm_analysis.getPathOutputTest("test_l1h_list.bin")

    # Copy settings file so that it is in the same place as the A matrix.
    settings_data = storm_analysis.getData("test/data/test_l1h.xml")
    settings_output = storm_analysis.getPathOutputTest("test_l1h.xml")
    shutil.copyfile(settings_data, settings_output)

    storm_analysis.removeFile(hres)
    storm_analysis.removeFile(mlist)

    from storm_analysis.L1H.cs_analysis import analyze
    analyze(movie_name, settings_output, hres, mlist)
    

if (__name__ == "__main__"):
    test_l1h()

