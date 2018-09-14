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

def test_get_attr_1():
    """
    Test no such parameter.
    """
    # Load some parameters.
    original = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    p1 = params.ParametersDAO().initFromFile(original, warnings = True)

    try:
        v1 = p1.getAttr("foo")
    except params.ParametersException:
        return

    assert False

def test_get_attr_2():
    """
    Test no such parameter with default value.
    """
    # Load some parameters.
    original = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    p1 = params.ParametersDAO().initFromFile(original, warnings = True)

    try:
        v1 = p1.getAttr("foo", default = "bar")
    except params.ParametersException:
        return

    assert False    

def test_get_attr_3():
    """
    Test getting parameters.
    """
    # Load some parameters.
    original = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    p1 = params.ParametersDAO().initFromFile(original, warnings = True)

    v1 = p1.getAttr("start_frame")
    v1 = p1.getAttr("x_start", 1)

def test_get_help_1():
    """
    Test getting help.
    """
    # Load some parameters.
    original = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    p1 = params.ParametersDAO().initFromFile(original, warnings = True)

    v1 = p1.helpAttr("max_frame")
    v1 = p1.helpAttr("convert_to")

def test_get_help_2():
    """
    Test getting help with a parameter that does not exist.
    """
    # Load some parameters.
    original = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    p1 = params.ParametersDAO().initFromFile(original, warnings = True)

    try:
        v1 = p1.helpAttr("foo")
    except params.ParametersException:
        return

    assert False

def test_pretty_print_1():
    # This just tests that it doesn't fail. It does not check the
    # formatting of the output file.

    original = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    output = storm_analysis.getPathOutputTest("test_pp1.xml")

    # Load parameters.
    p1 = params.ParametersDAO().initFromFile(original, warnings = True)

    # Convert back to XML.
    p1.toXMLFile(output, pretty = True)

def test_change_attr_1():
    """
    Test not specifying the node type.
    """
    p1 = params.ParametersDAO()

    try:
        v1 = p1.changeAttr("find_max_radius", 1)
    except params.ParametersException:
        return

    assert False

def test_change_attr_2():
    """
    Test specifying a node type that doesn't exist.
    """
    p1 = params.ParametersDAO()
        
    try:
        v1 = p1.changeAttr("find_max_radius", 1, node_type = "foo")
    except params.ParametersException:
        return

    assert False    
    
    
if (__name__ == "__main__"):
    test_xml_round_trip()
    test_get_attr_1()
    test_get_attr_2()
    test_get_attr_3()
    test_get_help_1()
    test_get_help_2()
    test_pretty_print_1()
    test_change_attr_1()
    test_change_attr_2()
