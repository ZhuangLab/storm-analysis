#!/usr/bin/env python
"""
Test Insight3 file format IO.
"""
import numpy

from xml.etree import ElementTree

import storm_analysis

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

def test_i3dtype_1():
    """
    Test conversion to and from the fitter format.
    """
    x_size = 100
    y_size = 100
    frame = 10
    nm_per_pixel = 100.0
    
    data_in = i3dtype.createDefaultI3Data(10)
    i3dtype.posSet(data_in, 'x', numpy.arange(10))
    i3dtype.posSet(data_in, 'y', numpy.arange(10) + 30.0)
    i3dtype.posSet(data_in, 'z', numpy.arange(10) + 60.0)
    i3dtype.setI3Field(data_in, 'fr', frame)

    peaks = i3dtype.convertToMultiFit(data_in, frame, nm_per_pixel)
    data_out = i3dtype.createFromMultiFit(peaks, frame, nm_per_pixel)

    fields = ['x', 'ax', 'w']
    for i in range(10):
        for field in fields:
            assert(abs(data_in[field][i] - data_out[field][i]) < 1.0e-6)
    
def test_i3dtype_2():
    """
    Test creation and setting position with an array.
    """
    i3_locs = i3dtype.createDefaultI3Data(10)
    i3dtype.posSet(i3_locs, 'x', numpy.arange(10))

    for i in range(10):
        assert(abs(i3_locs["x"][i] - i) < 1.0e-6)
        assert(abs(i3_locs["xc"][i] - i) < 1.0e-6)

def test_i3dtype_3():
    """
    Test creation and setting with a scalar.
    """
    i3_locs = i3dtype.createDefaultI3Data(10)
    i3dtype.posSet(i3_locs, 'x', 10.0)

    for i in range(10):
        assert(abs(i3_locs["x"][i] - 10.0) < 1.0e-6)
        assert(abs(i3_locs["xc"][i] - 10.0) < 1.0e-6)

def test_write_read_1():
    """
    Test writing and reading.
    """
    bin_name = storm_analysis.getPathOutputTest("test_insight3io.bin")

    i3_locs = i3dtype.createDefaultI3Data(10)
    i3dtype.posSet(i3_locs, 'x', 10.0)

    with writeinsight3.I3Writer(bin_name) as i3:
        i3.addMolecules(i3_locs)

    i3_in = readinsight3.loadI3File(bin_name, verbose = False)
    assert(numpy.allclose(i3_locs['x'], i3_in['x']))

def test_write_read_2():
    """
    Test I3Reader.
    """
    bin_name = storm_analysis.getPathOutputTest("test_insight3io.bin")

    i3_locs = i3dtype.createDefaultI3Data(10)
    i3dtype.posSet(i3_locs, 'x', 10.0)

    with writeinsight3.I3Writer(bin_name) as i3:
        i3.addMolecules(i3_locs)

    i3_reader = readinsight3.I3Reader(bin_name)

    # Read localizations.
    i3_in = i3_reader.nextBlock()
    assert(numpy.allclose(i3_locs['x'], i3_in['x']))
    assert(i3_in is not False)

    # This should return False because there are no more localizations.
    i3_in = i3_reader.nextBlock()
    assert(i3_in is False)

def test_write_read_3():
    """
    Test I3Reader on an empty file.
    """
    bin_name = storm_analysis.getPathOutputTest("test_insight3io.bin")

    i3_locs = i3dtype.createDefaultI3Data(10)
    i3dtype.posSet(i3_locs, 'x', 10.0)

    with writeinsight3.I3Writer(bin_name) as i3:
        pass

    i3_reader = readinsight3.I3Reader(bin_name)

    # Read localizations.
    i3_in = i3_reader.nextBlock()
    assert(i3_in is False)
    
def test_good_i3():
    mlist_name = storm_analysis.getPathOutputTest("test_i3_io_mlist.bin")

    # Create data.
    locs = i3dtype.createDefaultI3Data(100)

    # Save the data.
    with writeinsight3.I3Writer(mlist_name) as i3w:
        i3w.addMolecules(locs)

    # Read the data.
    locs = readinsight3.loadI3File(mlist_name)
    assert(locs.shape[0] == 100)

def test_good_i3_metadata():
    mlist_name = storm_analysis.getPathOutputTest("test_i3_io_mlist.bin")

    # Create data.
    locs = i3dtype.createDefaultI3Data(100)

    # Save data and metadata.
    i3w = writeinsight3.I3Writer(mlist_name)
    i3w.addMolecules(locs)
    etree = ElementTree.Element("xml")
    test = ElementTree.SubElement(etree, "test")
    test.text = "test"
    i3w.closeWithMetadata(ElementTree.tostring(etree, 'ISO-8859-1'))

    # Read the data.
    locs = readinsight3.loadI3File(mlist_name)
    assert(locs.shape[0] == 100)

    # Read the metadata.
    metadata = readinsight3.loadI3Metadata(mlist_name)
    assert(metadata.find("test").text == "test")

def test_bad_i3():
    mlist_name = storm_analysis.getPathOutputTest("test_i3_io_mlist.bin")

    # Create data.
    locs = i3dtype.createDefaultI3Data(100)

    # Save the data.
    i3w = writeinsight3.I3Writer(mlist_name)
    i3w.addMolecules(locs)
    i3w.fp.close()

    # Read the data.
    locs = readinsight3.loadI3File(mlist_name)
    assert(locs is None)
    
if (__name__ == "__main__"):
    test_i3dtype_1()
    test_i3dtype_2()
    test_i3dtype_3()
    test_write_read_1()
    test_write_read_2()
    test_write_read_3()
    test_good_i3()
    test_good_i3_metadata()
    test_bad_i3()    
