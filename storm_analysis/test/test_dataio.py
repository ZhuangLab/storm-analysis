#!/usr/bin/env python
"""
Test SMLM movie IO.
"""
import numpy
import tifffile

import storm_analysis

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.datawriter as datawriter


def test_io_1():
    """
    Test DAX movie IO.
    """
    movie_h = 50
    movie_w = 40
    movie_l = 10
    
    data = numpy.random.randint(0, 60000, (movie_h, movie_w)).astype(numpy.uint16)

    movie_name = storm_analysis.getPathOutputTest("test_dataio.dax")

    # Write dax movie.
    wr = datawriter.inferWriter(movie_name)
    for i in range(movie_l):
        wr.addFrame(data)
    wr.close()
        
    # Read & check.
    rd = datareader.inferReader(movie_name)
    [mw, mh, ml] = rd.filmSize()
    
    assert(mh == movie_h)
    assert(mw == movie_w)
    assert(ml == movie_l)
    assert(numpy.allclose(data, rd.loadAFrame(0)))

    rd.close()

def test_io_2():
    """
    Test TIF movie IO.
    """
    movie_h = 50
    movie_w = 40
    movie_l = 10
    
    data = numpy.random.randint(0, 60000, (movie_h, movie_w)).astype(numpy.uint16)

    movie_name = storm_analysis.getPathOutputTest("test_dataio.tif")

    # Write tif movie.
    wr = datawriter.inferWriter(movie_name)
    for i in range(movie_l):
        wr.addFrame(data)
    wr.close()
        
    # Read & check.
    rd = datareader.inferReader(movie_name)
    [mw, mh, ml] = rd.filmSize()

    assert(mh == movie_h)
    assert(mw == movie_w)
    assert(ml == movie_l)
    assert(numpy.allclose(data, rd.loadAFrame(0)))

def test_io_3():
    """
    Test FITS movie IO.
    """
    movie_h = 50
    movie_w = 40
    movie_l = 10
    
    data = numpy.random.randint(0, 60000, (movie_h, movie_w)).astype(numpy.uint16)

    movie_name = storm_analysis.getPathOutputTest("test_dataio.fits")

    # Write FITS movie.
    wr = datawriter.inferWriter(movie_name)
    for i in range(movie_l):
        wr.addFrame(data)
    wr.close()
        
    # Read & check.
    rd = datareader.inferReader(movie_name)
    [mw, mh, ml] = rd.filmSize()

    assert(mh == movie_h)
    assert(mw == movie_w)
    assert(ml == movie_l)
    assert(numpy.allclose(data, rd.loadAFrame(0)))
    
def test_io_4():
    """
    Test TIF movie IO (1 page, multiple frames per page).
    """
    movie_h = 50
    movie_w = 40
    movie_l = 10
    
    data = numpy.random.randint(0, 60000, (movie_l, movie_h, movie_w)).astype(numpy.uint16)

    movie_name = storm_analysis.getPathOutputTest("test_dataio.tif")

    # Write tif movie.
    with tifffile.TiffWriter(movie_name, imagej = True) as tf:
        tf.write(data, truncate = True)

    # Read & check.
    rd = datareader.inferReader(movie_name)
    [mw, mh, ml] = rd.filmSize()

    assert(mh == movie_h)
    assert(mw == movie_w)
    assert(ml == movie_l)
    for i in range(movie_l):
        assert(numpy.allclose(data[i,:,:], rd.loadAFrame(i)))

def test_io_5():
    """
    Test TIF movie IO (1 frame, 1 page).
    """
    movie_h = 50
    movie_w = 40
    movie_l = 1

    movie_name = storm_analysis.getPathOutputTest("test_dataio.tif")
    
    ## Standard Tiff.
    data = numpy.random.randint(0, 60000, (movie_h, movie_w)).astype(numpy.uint16)

    # Write tif movie.
    wr = datawriter.inferWriter(movie_name)
    for i in range(movie_l):
        wr.addFrame(data)
    wr.close()
        
    # Read & check.
    rd = datareader.inferReader(movie_name)
    [mw, mh, ml] = rd.filmSize()

    assert(mh == movie_h)
    assert(mw == movie_w)
    assert(ml == movie_l)
    assert(numpy.allclose(data, rd.loadAFrame(0)))

    ## 'imagej' Tiff.
    data = numpy.random.randint(0, 60000, (movie_l, movie_h, movie_w)).astype(numpy.uint16)

    movie_name = storm_analysis.getPathOutputTest("test_dataio.tif")

    # Write tif movie.
    with tifffile.TiffWriter(movie_name, imagej = True) as tf:
        tf.write(data, truncate = True)

    # Read & check.
    rd = datareader.inferReader(movie_name)
    [mw, mh, ml] = rd.filmSize()

    assert(mh == movie_h)
    assert(mw == movie_w)
    assert(ml == movie_l)
    assert(numpy.allclose(data[0,:,:], rd.loadAFrame(0)))

    
def test_io_6():
    """
    Test the dummyDaxFile() and singleFrameDax() helpers on a non-square
    movie.

    dummyDaxFile() built its frame as [x_size, y_size] while telling the
    writer that x_size was the width, so addFrame() asserted on anything
    that was not square. singleFrameDax(), immediately below it, has always
    had this right.
    """
    movie_w = 16
    movie_h = 8

    ## dummyDaxFile.
    movie_name = storm_analysis.getPathOutputTest("test_dataio_dummy.dax")

    datawriter.dummyDaxFile(movie_name, movie_w, movie_h)

    rd = datareader.inferReader(movie_name)
    [mw, mh, ml] = rd.filmSize()

    assert(mw == movie_w)
    assert(mh == movie_h)
    assert(ml == 1)

    frame = rd.loadAFrame(0)
    assert(frame.shape == (movie_h, movie_w))
    assert(numpy.allclose(frame, numpy.ones((movie_h, movie_w))))

    ## singleFrameDax, the same convention from the other side.
    data = numpy.random.randint(0, 60000, (movie_h, movie_w)).astype(numpy.uint16)

    movie_name = storm_analysis.getPathOutputTest("test_dataio_single.dax")

    datawriter.singleFrameDax(movie_name, data)

    rd = datareader.inferReader(movie_name)
    [mw, mh, ml] = rd.filmSize()

    assert(mw == movie_w)
    assert(mh == movie_h)
    assert(ml == 1)
    assert(numpy.allclose(data, rd.loadAFrame(0)))

def test_io_7():
    """
    DaxReader on an incomplete .inf file, and on a missing .dax.

    'number of frames' and the endianness had no fallback, so a .inf file
    without them gave AttributeError from inside loadAFrame(). A missing
    .dax constructed a reader that then failed with AttributeError on a
    None file pointer.
    """
    import pytest

    movie_h = 8
    movie_w = 16
    movie_l = 5

    data = numpy.random.randint(0, 60000, (movie_h, movie_w)).astype(numpy.uint16)

    movie_name = storm_analysis.getPathOutputTest("test_dataio_inf.dax")
    inf_name = storm_analysis.getPathOutputTest("test_dataio_inf.inf")

    wr = datawriter.inferWriter(movie_name)
    for i in range(movie_l):
        wr.addFrame(data)
    wr.close()

    with open(inf_name) as fp:
        inf_lines = fp.readlines()

    def writeInf(drop):
        with open(inf_name, "w") as fp:
            for line in inf_lines:
                if (drop is None) or (not drop in line):
                    fp.write(line)

    ## No 'number of frames', the file size gives it instead.
    writeInf("number of frames")
    rd = datareader.inferReader(movie_name)
    assert(rd.filmSize() == [movie_w, movie_h, movie_l])
    assert(numpy.allclose(data, rd.loadAFrame(0)))
    rd.close()

    ## No endianness, little endian is assumed.
    writeInf("data type")
    rd = datareader.inferReader(movie_name)
    assert(numpy.allclose(data, rd.loadAFrame(0)))
    rd.close()

    ## A complete .inf but no .dax, and the error names the file.
    missing_name = storm_analysis.getPathOutputTest("test_dataio_missing.dax")
    storm_analysis.removeFile(missing_name)
    with open(storm_analysis.getPathOutputTest("test_dataio_missing.inf"), "w") as fp:
        for line in inf_lines:
            fp.write(line)

    with pytest.raises(IOError) as err:
        datareader.inferReader(missing_name)
    assert("test_dataio_missing.dax" in str(err.value))


if (__name__ == "__main__"):
    test_io_1()
    test_io_2()
    test_io_3()
    test_io_6()
    test_io_7()
    
