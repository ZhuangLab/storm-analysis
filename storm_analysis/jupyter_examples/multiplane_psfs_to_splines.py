#!/usr/bin/env python
"""
Helper functions for Multiplane PSF to spline conversion.

Hazen 10/17
"""
import storm_analysis.sa_library.parameters as parameters

pixel_size = 100.0
spline_z_range = 0.75
z_value = [-0.3, 0.0, 0.3]

def multiplaneXML():
    """
    Create a Multiplane parameters object.
    """
    params = parameters.ParametersMultiplaneArb()

    params.setAttr("max_frame", "int", -1)    
    params.setAttr("start_frame", "int", -1)
    
    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("find_max_radius", "int", 2)
    params.setAttr("independent_heights", "int", 0)
    params.setAttr("iterations", "int", 20)
    params.setAttr("mapping", "filename", "map.map")
    params.setAttr("no_fitting", "int", 0)
    params.setAttr("pixel_size", "float", pixel_size)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", 6.0)
    params.setAttr("weights", "filename", "weights.npy")
    params.setAttr("z_value", "float-array", z_value)

    params.setAttr("channel0_cal", "filename", "calib.npy")
    params.setAttr("channel1_cal", "filename", "calib.npy")

    params.setAttr("channel0_ext", "string", "_c1.dax")
    params.setAttr("channel1_ext", "string", "_c2.dax")

    params.setAttr("channel0_offset", "int", 0)
    params.setAttr("channel1_offset", "int", 0)

    params.setAttr("spline0", "filename", "c1_psf.spline")
    params.setAttr("spline1", "filename", "c2_psf.spline")

    # Don't do tracking.
    params.setAttr("descriptor", "string", "1")
    params.setAttr("radius", "float", "0.0")

    params.setAttr("max_z", "float", str(spline_z_range + 0.001))
    params.setAttr("min_z", "float", str(-(spline_z_range - 0.001)))

    # Don't do drift-correction.
    params.setAttr("d_scale", "int", 2)
    params.setAttr("drift_correction", "int", 0)
    params.setAttr("frame_step", "int", 500)
    params.setAttr("z_correction", "int", 0)

    params.toXMLFile("multiplane.xml")


if (__name__ == "__main__"):
    multiplaneXML()
