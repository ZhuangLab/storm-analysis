#!/usr/bin/env python

import storm_analysis

import storm_analysis.sa_library.parameters as params

def test_xml_round_trip():

    original = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    output = storm_analysis.getPathOutputTest("test_rt.xml")

    # Load parameters.
    p1 = params.ParametersDAO().initFromFile(original, warnings = True)

    # Convert back to XML.
    p1.toXMLFile(output)

    # And read again.
    p2 = params.ParametersDAO().initFromFile(output, warnings = True)

    # And compare.
    if (p1.getAttr("threshold") != p2.getAttr("threshold")):
        raise Exception("Parameter XML input / output failed.")


if (__name__ == "__main__"):
    test_xml_round_trip()
