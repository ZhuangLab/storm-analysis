#!/usr/bin/env python
"""
Test of sa_library.analysis_io
"""
import numpy

import storm_analysis
import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sCMOS.reslice_calibration as resliceCalibration


cal_size = (20, 18)


def test_cal_v0():
    """
    Test loading a v0 calibration file.
    """
    cal_file = storm_analysis.getPathOutputTest("calib.npy")

    offset = numpy.random.uniform(size = cal_size)
    variance = numpy.random.uniform(size = cal_size)
    gain = numpy.random.uniform(size = cal_size)

    numpy.save(cal_file,
               [numpy.transpose(offset),
                numpy.transpose(variance),
                numpy.transpose(gain)])

    [o, v, g, r] = analysisIO.loadCMOSCalibration(cal_file)
    assert(numpy.allclose(offset, o))
    assert(numpy.allclose(variance, v))
    assert(numpy.allclose(gain, g))
    assert(numpy.allclose(r, numpy.ones(cal_size)))


def test_cal_v1():
    """
    Test loading a v1 calibration file.
    """
    cal_file = storm_analysis.getPathOutputTest("calib.npy")

    offset = numpy.random.uniform(size = cal_size)
    variance = numpy.random.uniform(size = cal_size)
    gain = numpy.random.uniform(size = cal_size)
    rqe = numpy.random.uniform(size = cal_size)

    numpy.save(cal_file, [offset, variance, gain, 1])

    [o, v, g, r] = analysisIO.loadCMOSCalibration(cal_file)
    assert(numpy.allclose(offset, o))
    assert(numpy.allclose(variance, v))
    assert(numpy.allclose(gain, g))
    assert(numpy.allclose(r, numpy.ones(cal_size)))
        

def test_cal_v2():
    """
    Test loading a v2 calibration file.
    """
    cal_file = storm_analysis.getPathOutputTest("calib.npy")

    offset = numpy.random.uniform(size = cal_size)
    variance = numpy.random.uniform(size = cal_size)
    gain = numpy.random.uniform(size = cal_size)
    rqe = numpy.random.uniform(size = cal_size)

    numpy.save(cal_file, [offset, variance, gain, rqe, 2])

    [o, v, g, r] = analysisIO.loadCMOSCalibration(cal_file)
    assert(numpy.allclose(offset, o))
    assert(numpy.allclose(variance, v))
    assert(numpy.allclose(gain, g))
    assert(numpy.allclose(rqe, r))

    
def test_cal_error_handling():
    """
    Test calibration file loader error handling.
    """
    cal_file = storm_analysis.getPathOutputTest("calib.npy")

    offset = numpy.random.uniform(size = cal_size)

    okay = False
    numpy.save(cal_file, [offset, offset, offset, 2])
    try:
        [o, v, g, r] = analysisIO.loadCMOSCalibration(cal_file)
    except analysisIO.AnalysisIOException:
        okay = True
    assert okay

    okay = False
    numpy.save(cal_file, [offset, offset, offset, offset, 1])
    try:
        [o, v, g, r] = analysisIO.loadCMOSCalibration(cal_file)
    except analysisIO.AnalysisIOException:
        okay = True
    assert okay


def test_cal_reslice():
    """
    Test that calibration re-slicing works as expected.
    """
    cal_file = storm_analysis.getPathOutputTest("calib.npy")
    cal_rs_file = storm_analysis.getPathOutputTest("calib_rs.npy")

    offset = numpy.random.uniform(size = cal_size)
    variance = numpy.random.uniform(size = cal_size)
    gain = numpy.random.uniform(size = cal_size)
    rqe = numpy.random.uniform(size = cal_size)

    numpy.save(cal_file, [offset, variance, gain, rqe, 2])
    resliceCalibration.resliceCalibration(cal_file, cal_rs_file, 1, 2, 10, 11)
    
    [o, v, g, r] = analysisIO.loadCMOSCalibration(cal_rs_file)
    assert(o.shape == (11,10))
    assert(numpy.allclose(offset[2:13,1:11], o))
    assert(numpy.allclose(variance[2:13,1:11], v))
    assert(numpy.allclose(gain[2:13,1:11], g))
    assert(numpy.allclose(rqe[2:13,1:11], r))

    assert(analysisIO.isLatestFormat(cal_rs_file))
    
    
if (__name__ == "__main__"):
    test_cal_v0()
    test_cal_v1()
    test_cal_v2()
    test_cal_error_handling()
    test_cal_reslice()
