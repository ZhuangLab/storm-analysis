#!/usr/bin/env python
"""
Test of multi_plane.analysis_io
"""

import numpy
import os
import tifffile

import storm_analysis
import storm_analysis.sa_library.parameters as params

import storm_analysis.multi_plane.analysis_io as analysis_io


im_size = (20, 18)


def configureTest():
    """
    These tests have a lot of setup. This function takes care of this.
    """
    mparams = params.ParametersMultiplane()
        
    # sCMOS calibration files.
    gain = numpy.ones(im_size)
    offset = numpy.zeros(im_size)
    variance = numpy.ones(im_size)
    rqe = numpy.ones(im_size)
    
    cal_file = storm_analysis.getPathOutputTest("c1_cal.npy")
    numpy.save(cal_file, [offset, variance, gain, rqe, 2])
    mparams.changeAttr("channel0_cal", cal_file)

    cal_file = storm_analysis.getPathOutputTest("c2_cal.npy")
    numpy.save(cal_file, [offset, variance, gain, rqe, 2])
    mparams.changeAttr("channel1_cal", cal_file)

    mparams.changeAttr("channel0_ext", "_c1.tif")
    mparams.changeAttr("channel1_ext", "_c2.tif")

    mparams.changeAttr("channel0_offset", 0)
    mparams.changeAttr("channel1_offset", 0)

    return mparams
    

def test_offset_1():
    """
    Movie reading test with zero offset.
    """
    # Create XML file and sCMOS calibration files.
    mparams = configureTest()
    xml_file = storm_analysis.getPathOutputTest("mp_analysis.xml")
    mparams.toXMLFile(xml_file)

    # Create movies.
    mv_file = storm_analysis.getPathOutputTest("mp_c1.tif")
    with tifffile.TiffWriter(mv_file) as tf:
        for i in range(10):
            im = (i+1) * numpy.ones(im_size, dtype = numpy.int16)
            tf.save(im.astype(numpy.int16))

    mv_file = storm_analysis.getPathOutputTest("mp_c2.tif")
    with tifffile.TiffWriter(mv_file) as tf:
        for i in range(10):
            im = (i+1) * numpy.ones(im_size, dtype = numpy.int16)
            tf.save(im.astype(numpy.int16))

    # Create and test MPMovieReader    
    mpmr = analysis_io.MPMovieReader(base_name = os.path.join(storm_analysis.getPathOutputTest(), "mp"),
                                     parameters = params.ParametersMultiplane().initFromFile(xml_file))
    mpmr.setup(0)
    
    [mx, my, ml] = mpmr.getFilmSize()
    assert(mx == im_size[1])
    assert(my == im_size[0])
    assert(ml == 10)

    # Check that all of the frames are the same.
    for i in range(ml):
        frames = mpmr.getFrames(i)
        assert(numpy.allclose(frames[0], frames[1]))
        

def test_offset_2():
    """
    Movie reading test with non-zero offset.
    """
    # Create XML file and sCMOS calibration files.
    mparams = configureTest()
    xml_file = storm_analysis.getPathOutputTest("mp_analysis.xml")

    mparams.changeAttr("channel1_offset", 1)
    mparams.toXMLFile(xml_file)

    # Create movies.
    mv_file = storm_analysis.getPathOutputTest("mp_c1.tif")
    with tifffile.TiffWriter(mv_file) as tf:
        for i in range(10):
            im = (i+1) * numpy.ones(im_size, dtype = numpy.int16)
            tf.save(im.astype(numpy.int16))

    mv_file = storm_analysis.getPathOutputTest("mp_c2.tif")
    with tifffile.TiffWriter(mv_file) as tf:
        for i in range(10):
            im = (i+1) * numpy.ones(im_size, dtype = numpy.int16)
            tf.save(im.astype(numpy.int16))

    # Create and test MPMovieReader    
    mpmr = analysis_io.MPMovieReader(base_name = os.path.join(storm_analysis.getPathOutputTest(), "mp"),
                                     parameters = params.ParametersMultiplane().initFromFile(xml_file))
    mpmr.setup(0)
    
    [mx, my, ml] = mpmr.getFilmSize()
    assert(mx == im_size[1])
    assert(my == im_size[0])
    assert(ml == 9)

    # Check that all of the frames are the same.
    for i in range(ml):
        frames = mpmr.getFrames(i)
        assert(numpy.allclose(frames[0], frames[1] - 1))

        
def test_offset_3():
    """
    Movie reading test with non-zero offset for channel 0.
    """
    # Create XML file and sCMOS calibration files.
    mparams = configureTest()
    xml_file = storm_analysis.getPathOutputTest("mp_analysis.xml")

    mparams.changeAttr("channel0_offset", 1)
    mparams.toXMLFile(xml_file)

    # Create movies.
    mv_file = storm_analysis.getPathOutputTest("mp_c1.tif")
    with tifffile.TiffWriter(mv_file) as tf:
        for i in range(10):
            im = (i+1) * numpy.ones(im_size, dtype = numpy.int16)
            tf.save(im.astype(numpy.int16))

    mv_file = storm_analysis.getPathOutputTest("mp_c2.tif")
    with tifffile.TiffWriter(mv_file) as tf:
        for i in range(10):
            im = (i+1) * numpy.ones(im_size, dtype = numpy.int16)
            tf.save(im.astype(numpy.int16))

    # Create and test MPMovieReader    
    mpmr = analysis_io.MPMovieReader(base_name = os.path.join(storm_analysis.getPathOutputTest(), "mp"),
                                     parameters = params.ParametersMultiplane().initFromFile(xml_file))
    mpmr.setup(0)
    
    [mx, my, ml] = mpmr.getFilmSize()
    assert(mx == im_size[1])
    assert(my == im_size[0])
    assert(ml == 9)

    # Check that all of the frames are the same.
    for i in range(ml):
        frames = mpmr.getFrames(i)
        assert(numpy.allclose(frames[0]-1, frames[1]))


if (__name__ == "__main__"):
    test_offset_1()
    test_offset_2()
    test_offset_3()


    
    
