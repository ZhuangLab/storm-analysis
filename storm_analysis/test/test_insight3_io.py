#!/usr/bin/env python

from xml.etree import ElementTree

import storm_analysis

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3


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
    test_good_i3()
    test_good_i3_metadata()
    test_bad_i3()
    
