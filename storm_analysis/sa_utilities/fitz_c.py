#!/usr/bin/env python
"""
Uses the fitz C library to calculate localization Z position
based on it's width in x and y. This is used by the '3d' model
of 3D-DAOSTORM / sCMOS analysis.

Note that because the widths in x/y are transposed in the HDF5
format relative to the Insight3 bin format you may need to 
update your calibration parameters.

Hazen 1/18
"""
import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os

import storm_analysis.sa_library.loadclib as loadclib
import storm_analysis.sa_library.sa_h5py as saH5Py

c_fitz = loadclib.loadCLibrary("storm_analysis.sa_utilities", "fitz")

c_fitz.cleanup.argtypes = [ctypes.c_void_p]

c_fitz.initialize.argtypes = [ndpointer(dtype=numpy.float64),
                              ndpointer(dtype=numpy.float64),
                              ctypes.c_double,
                              ctypes.c_double,
                              ctypes.c_double,
                              ctypes.c_double]
c_fitz.initialize.restype = ctypes.c_void_p

c_fitz.findBestZ.argtypes = [ctypes.c_void_p,
                             ctypes.c_double,
                             ctypes.c_double]
c_fitz.findBestZ.restype = ctypes.c_double


def calcSxSy(wx_params, wy_params, z):
    """
    Return sigma x and sigma y given the z calibration parameters.
    """
    zx = (z - wx_params[1])/wx_params[2]
    sx = 0.5 * wx_params[0] * numpy.sqrt(1.0 + zx*zx + wx_params[3]*zx*zx*zx + wx_params[4]*zx*zx*zx*zx)
    zy = (z - wy_params[1])/wy_params[2]
    sy = 0.5 * wy_params[0] * numpy.sqrt(1.0 + zy*zy + wy_params[3]*zy*zy*zy + wy_params[4]*zy*zy*zy*zy)
    return [sx, sy]


def fitz(h5_name, cutoff, pixel_size, wx_params, wy_params, z_min, z_max, z_step = 0.001):
    """
    This processes both the raw and the tracked localizations.

    cutoff - Max allowed distance from the wx/wy versus Z curve, units unclear.
    pixel_size - nanometers/pixel.
    wx_params, wy_params - These are in nanometers / dimensionless, as for
                           example the values returned by zee-calibrator in the
                           storm-control project.
    z_min, z_max - Minimum and maximum values in microns.
    z_step - Step size of Z search in microns.
    """
    # Fit raw localizations.
    fitzRaw(h5_name, cutoff, pixel_size, wx_params, wy_params, z_min, z_max, z_step)

    # Fit tracks.
    fitzTracks(h5_name, cutoff, pixel_size, wx_params, wy_params, z_min, z_max, z_step)
    

def fitzRaw(h5_name, cutoff, wx_params, wy_params, z_min, z_max, z_step):
    """
    This processes the raw localizations.

    Note: Localizations whose wx/wy values are too far from the calibration
          curve will be given a z value that is less than z_min.
    """
    zfit_data = c_fitz.initialize(numpy.ascontiguousarray(wx_params),
                                  numpy.ascontiguousarray(wy_params),
                                  z_min * 1000.0,
                                  z_max * 1000.0,
                                  z_step * 1000.0,
                                  cutoff)

    # Fit raw localizations & save z value (in microns).
    with saH5Py.SAH5Py(h5_name) as h5:
        pixel_size = h5.getPixelSize()
        for fnum, locs in h5.localizationsIterator():
            z_vals = numpy.zeros(locs["xsigma"].size, dtype = numpy.float64)
            for i in range(locs["xsigma"].size):
                wx = pixel_size * 2.0 * locs["xsigma"][i]
                wy = pixel_size * 2.0 * locs["ysigma"][i]
                z_vals[i] = c_fitz.findBestZ(zfit_data, wx, wy) * 1.0e-3
            h5.addLocalizationZ(z_vals, fnum)

    c_fitz.cleanup(zfit_data)


def fitzTracks(h5_name, cutoff, wx_params, wy_params, z_min, z_max, z_step):
    """
    This processes the tracked localizations.

    Note: Localizations whose wx/wy values are too far from the calibration
          curve will be given a z value that is less than z_min and also
          assigned to category 9.
    """
    zfit_data = c_fitz.initialize(numpy.ascontiguousarray(wx_params),
                                  numpy.ascontiguousarray(wy_params),
                                  z_min * 1000.0,
                                  z_max * 1000.0,
                                  z_step * 1000.0,
                                  cutoff)

    # Fit tracked localizations & save z value (in microns).
    with saH5Py.SAH5Py(h5_name) as h5:
        pixel_size = h5.getPixelSize()
        for index, locs in enumerate(h5.tracksIterator()):
            z_vals = numpy.zeros(locs["xsigma"].size, dtype = numpy.float64)
            for i in range(locs["xsigma"].size):
                wx = pixel_size * 2.0 * locs["xsigma"][i]/locs["track_length"][i]                    
                wy = pixel_size * 2.0 * locs["ysigma"][i]/locs["track_length"][i]
                z_vals[i] = c_fitz.findBestZ(zfit_data, wx, wy) * 1.0e-3
            h5.addTrackData(z_vals, index, "z")

    c_fitz.cleanup(zfit_data)

    

if (__name__ == "__main__"):
    
    import argparse

    import storm_analysis.sa_library.parameters as params

    parser = argparse.ArgumentParser(description = 'Z fitting given Wx, Wy calibration curves.')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations file.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()

    parameters = params.ParametersDAO().initFromFile(args.settings)

    [wx_params, wy_params] = parameters.getWidthParams()
    [min_z, max_z] = parameters.getZRange()
        
    fitz(args.mlist,
         parameters.getAttr("cutoff"),
         wx_params,
         wy_params,
         min_z,
         max_z,
         parameters.getAttr("z_step", 0.001))
