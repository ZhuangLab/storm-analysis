#!/usr/bin/env python
"""
Calculate z fitting values for 3D astigmatism imaging and the 
'3d' or 'Z' fitting model in 3D-DAOSTORM.

Internal calculations are done in pixels and microns. However 
in order to work correctly with the current analysis pipeline
they must be converted to nanometers.

Hazen 01/18
"""
import copy
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import scipy
import scipy.optimize

from xml.dom import minidom
from xml.etree import ElementTree

import storm_analysis
import storm_analysis.sa_library.sa_h5py as saH5Py


class ZCalibrationException(Exception):
    pass


def calibrate(hdf5, zoffsets, fit_order, outliers, no_plots = False, wx_wy_fields = ["xsigma", "ysigma"], stage_tilt = [0.0, 0.0, 0.0]):
    """
    Run all the steps in Z calibration in a single step.

    hdf5 - The HDF5 file.
    zoffsets - Text file containing z offsets.
    fit_order - An integer in the range 0-4.
    outliers - Sigma threshold for removing outliers in Wx, Wy.
    no_plots - Don't make a plot of the fit.
    wx_wy_fields - Which fields to use for wx and wy.
    stage_tilt - Array to use for correction of stage tilt (if any).
    """
        
    # Load the data.
    [wx, wy, z, pixel_size] = loadWxWyZData(hdf5, zoffsets, wx_wy_fields = wx_wy_fields, stage_tilt = stage_tilt)

    # Fit curves.
    print("Fitting (round 1).")
    [wx_params, wy_params] = fitDefocusingCurves(wx, wy, z, n_additional = fit_order)

    # Remove outliers.
    print("Removing outliers.")
    [t_wx, t_wy, t_z] = removeOutliers(wx, wy, z, wx_params, wy_params, outliers)

    # Redo fit.
    print("Fitting (round 2).")
    [wx_params, wy_params] = fitDefocusingCurves(t_wx, t_wy, t_z, n_additional = fit_order)

    # Plot fit.
    if not no_plots:
        plotFit(wx, wy, z, t_wx, t_wy, t_z, wx_params, wy_params)

    # Print results.
    prettyPrint(wx_params, wy_params, pixel_size)

    return [wx_params, wy_params]

    
def convertUnits(ww_params, pixel_size):
    """
    Convert to the units currently used by the analysis pipeline. This also
    adds additional zeros as needed to match the length expected by the
    analysis pipeline.
    """
    ww_params = copy.copy(ww_params)
    
    ww_params[0] = ww_params[0] * pixel_size
    ww_params[1] = ww_params[1] * 1.0e+3
    ww_params[2] = ww_params[2] * 1.0e+3

    for i in range(len(ww_params), 7):
        ww_params.append(0.0)
        
    return ww_params

    
def doFit(w, z, z_params, n_additional):
    """
    w - Numpy array containing localizations widths in pixels.
    z - Numpy array containing localization z positions in microns.
    z_params - Initial guess for the fitting parameters [pixels, microns, microns].
    n_additional - The number of additional fitting parameters to use,
                   these are 'A,B,C,D'.
    """
    assert(w.size == z.size)
    assert(len(z_params) == 3)
    assert(n_additional >= 0)
    assert(n_additional < 5)

    # Fit zero order to get initial estimate.
    zfit = doSingleFit(w, z, z_params, 0)

    # If this fails, try inverting the second term.
    if zfit is None:
        z_params[1] = -z_params[1]
        doSingleFit(w, z, z_params, 0)

    # Throw exception if we still failed.
    if zfit is None:
        raise ZCalibrationException("Defocus curve fitting failed (order 0).")

    # Fit again with requested number of additional parameters.
    z_params = zfit.tolist()
    zfit = doSingleFit(w, z, z_params, n_additional)

    # Throw exception if we failed.
    if zfit is None:
        raise ZCalibrationException("Defocus curve fitting failed (order 0).")

    return zfit.tolist()


def doSingleFit(w, z, z_params, n_additional):
    """
    w - Numpy array containing localizations widths in pixels.
    z - Numpy array containing localization z positions in microns.
    z_params - Initial guess for the fitting parameters [pixels, microns, microns].
    n_additional - The number of additional fitting parameters to use,
                   these are 'A,B,C,D'.
    """
    def zcalib_fitter(fit_p, w, z):
        return zcalib_fitters[n_additional](fit_p, z) - w

    # Add additional fitting parameters.
    for i in range(n_additional):
        z_params.append(0.0)
    [results, success] = scipy.optimize.leastsq(zcalib_fitter, z_params, args=(w, z))
    if (success < 1) or (success > 4):
        return None
    else:
        return results


def fitDefocusingCurves(wx, wy, z, n_additional = 0, z_params = None):
    """
    wx - Numpy array with widths in x in pixels.
    wy - Numpy array with widths in y in pixels.
    z - Numpy array with z values in microns.
    n_additional - The number of additional fitting parameters to use,
                   these are 'A,B,C,D'.
    z_params - (Optional) 3 element list containing an initial guess
               for the z fitting parameters.
    """
    if z_params is not None:
        assert(len(z_params) == 3)
    else:
        z_params = [3.0, 0.3, 0.5]

    # Fit for defocusing parameters.
    wx_params = doFit(wx, z, z_params, n_additional)
    wy_params = doFit(wy, z, z_params, n_additional)
    
    return [wx_params, wy_params]


def fitTilt(h5_name, start = 0, stop = 1):
    """
    h5_name - The name of the HDF5 localization file (must be fit for z).

    Use the localizations in the range start <= frame number < stop. This 
    should be a region of the movie where the stage is near zero and also
    not moving.
    """

    # Load localizations.
    with saH5Py.SAH5Py(h5_name) as h5:
        locs = h5.getLocalizationsInFrameRange(start, stop, fields = ["x", "y", "z"])
        
    # Find the best fit plane through x,y,z.
    def fitfn(p):
        zf = p[0] + p[1]*locs["x"] + p[2]*locs["y"]
        return locs["z"] - zf
    
    params = [scipy.mean(locs["z"]), 0.0, 0.0]
    [results, success] = scipy.optimize.leastsq(fitfn, params)
    
    if (success < 1) or (success > 4):
        raise ZCalibrationException("fitTilt: fit failed!")

    return results


def loadWxWyZData(h5_name, zfile_name, wx_wy_fields = ["xsigma", "ysigma"], stage_tilt = [0.0, 0.0, 0.0]):
    """
    h5_name - The name of the HDF5 localization file.
    zfile_name - The name of the text file containing z offset data.
    """
    wx_field = wx_wy_fields[0]
    wy_field = wx_wy_fields[1]
    
    # This file contains two columns, the first is whether or not
    # the data in this frame should be used (0 = No, 1 = Yes) and
    # the second contains the z offset in microns.
    #
    # For movies acquired using storm-control this can be created
    # from the .off file using storm_analysis/spliner/offset_to_z.py.
    #
    z_data = numpy.loadtxt(zfile_name, ndmin = 2)

    # Create arrays with wx, wy, z data.
    wx = None
    wy = None
    z = None
    with saH5Py.SAH5Py(h5_name) as h5:
        pixel_size = h5.getPixelSize()
        for curf, locs in h5.localizationsIterator(fields = wx_wy_fields + ["x", "y"]):
            if (int(z_data[curf,0]) == 0):
                continue

            # Calculate stage tilt corrected z values.
            tz = tiltCompensation(locs["x"],
                                  locs["y"],
                                  numpy.ones(locs["x"].size) * z_data[curf, 1],
                                  stage_tilt)

            if wx is None:
                wx = 2.0 * locs[wx_field]
                wy = 2.0 * locs[wy_field]
                z = tz
            else:
                wx = numpy.concatenate((wx, 2.0 * locs[wx_field]))
                wy = numpy.concatenate((wy, 2.0 * locs[wy_field]))
                z = numpy.concatenate((z, tz))

    return [wx, wy, z, pixel_size]


def plotFit(wx, wy, z, t_wx, t_wy, t_z, wx_params, wy_params, z_range = 0.6):
    """
    wx - Numpy array with widths in x in pixels.
    wy - Numpy array with widths in y in pixels.
    z - Numpy array with z values in microns.
    wx_params - Parameters for wx versus z fit.
    wy_params - Parameters for wy versus z fit.
    """
    storm_analysis.configureMatplotlib()

    # Add a little dither to the x position.
    pz = z + numpy.random.normal(scale = 0.003, size = z.size)
    t_pz = t_z + numpy.random.normal(scale = 0.003, size = t_z.size)
    
    z_fit = numpy.arange(-z_range, z_range + 0.001, 0.01)

    zfn = zcalib_fitters[len(wx_params) - 3]
    wx_fit = zfn(wx_params, z_fit)
    wy_fit = zfn(wy_params, z_fit)

    # Plot discarded points in grey.
    pyplot.scatter(pz, wx, color = 'lightgray', s = 1)
    pyplot.scatter(pz, wx, color = 'lightgray', s = 1)

    # Plot points and fit.
    p1 = pyplot.scatter(t_pz, t_wx, color = 'r', label = 'Wx', s = 1)
    p2 = pyplot.scatter(t_pz, t_wy, color = 'g', label = 'Wy', s = 1)
    pyplot.plot(z_fit, wx_fit, color = 'black')
    pyplot.plot(z_fit, wy_fit, color = 'black')
    legend = pyplot.legend(handles = [p1, p2], loc=1)
    legend.get_frame().set_linewidth(2)
    legend.get_frame().set_edgecolor('black')
    pyplot.xlabel("microns")
    pyplot.ylabel("pixels")
    pyplot.show()


def prettyPrint(wx_params, wy_params, pixel_size):
    """
    Print parameters in a storm-analysis XML friendly format.

    wx_params - Parameters for wx versus z fit.
    wy_params - Parameters for wy versus z fit.
    pixel_size - Pixel size in nanometers.
    """
    wx_params = convertUnits(wx_params, pixel_size)
    wy_params = convertUnits(wy_params, pixel_size)

    # Create XML.
    etree = ElementTree.Element("xml")
    for i, elt in enumerate([wx_params, wy_params]):
        txt1 = "wx"
        if (i == 1):
            txt1 = "wy"
        for j, txt2 in enumerate(["_wo", "_c", "_d", "A", "B", "C", "D"]):
            tmp = ElementTree.SubElement(etree, txt1 + txt2)
            tmp.text = "{0:0.3f}".format(elt[j])

    # Pretty print it.
    reparsed = minidom.parseString(ElementTree.tostring(etree))
    print(reparsed.toprettyxml(indent = "   ", encoding = "ISO-8859-1").decode())


def removeOutliers(wx, wy, z, wx_params, wy_params, threshold):
    """
    This removes all wx, wy that are more than threshold sigma from the fit curve.
    """
    zfn = zcalib_fitters[len(wx_params) - 3]
    wx_fit = zfn(wx_params, z)
    wy_fit = zfn(wy_params, z)

    wx_dev = numpy.std(wx_fit - wx)
    wy_dev = numpy.std(wy_fit - wy)

    wx_mask = (numpy.abs(wx_fit - wx) < (threshold * wx_dev))
    wy_mask = (numpy.abs(wy_fit - wy) < (threshold * wy_dev))

    mask = numpy.logical_and(wx_mask, wy_mask)

    return [wx[mask], wy[mask], z[mask]]


def setWxWyParams(params, wx_params, wy_params, pixel_size):
    """
    Set calibration parameters based on wx_params and wy_params.

    params - A parameters object.
    wx_params - Parameters for wx versus z fit.
    wy_params - Parameters for wy versus z fit.
    pixel_size - Pixel size in nanometers.
    """
    wx_params = convertUnits(wx_params, pixel_size)
    wy_params = convertUnits(wy_params, pixel_size)

    # Create XML.
    for i, elt in enumerate([wx_params, wy_params]):
        txt1 = "wx"
        if (i == 1):
            txt1 = "wy"
        for j, txt2 in enumerate(["_wo", "_c", "_d", "A", "B", "C", "D"]):
            params.setAttr(txt1 + txt2, "float", elt[j])


def tiltCompensation(x, y, z, stage_tilt):
    """
    Corrects z for stage tilt.

    This is pretty simple, but it is factored out to make for the
    benefit of external modules.
    """
    return z + stage_tilt[0] + stage_tilt[1]*x + stage_tilt[2]*y

                           
# These are the different defocus curve fitting functions.
#
def zcalib0(p, z):
    """
    Z calibration fitting function with no additional parameters.
    """
    wo,c,d = p
    X = (z-c)/d
    return wo*numpy.sqrt(1.0 + numpy.power(X,2))

def zcalib1(p, z):
    """
    Z calibration fitting function with 1 additional parameters.
    """
    wo,c,d,A = p
    X = (z-c)/d
    return wo*numpy.sqrt(1.0 + numpy.power(X,2) + A * numpy.power(X,3))

def zcalib2(p, z):
    """
    Z calibration fitting function with 2 additional parameters.
    """
    wo,c,d,A,B = p
    X = (z-c)/d
    return wo*numpy.sqrt(1.0 + numpy.power(X,2) + A * numpy.power(X,3) + B * numpy.power(X,4))

def zcalib3(p, z):
    """
    Z calibration fitting function with 3 additional parameters.
    """
    wo,c,d,A,B,C = p
    X = (z-c)/d
    return wo*numpy.sqrt(1.0 + numpy.power(X,2) + A * numpy.power(X,3) + B * numpy.power(X,4) + C * numpy.power(X,5))

def zcalib4(p, z):
    """
    Z calibration fitting function with 4 additional parameters.
    """
    wo,c,d,A,B,C,D = p
    X = (z-c)/d
    return wo*numpy.sqrt(1.0 + numpy.power(X,2) + A * numpy.power(X,3) + B * numpy.power(X,4) + C * numpy.power(X,5) + D * numpy.power(X,6))

zcalib_fitters = [zcalib0, zcalib1, zcalib2, zcalib3, zcalib4]


if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = 'Fit for defocusing curve parameters.')

    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the localizations HDF5 file.")
    parser.add_argument('--zoffsets', dest='zoffsets', type=str, required=True,
                        help = "A text file with two space separated numbers on each line, the first is 1 of the frame is valid, 0 otherwise and the second is the z offset of the frame relative to the focal plane in microns.")
    parser.add_argument('--fit_order', dest='fit_order', type=int, required=False, default=2,
                        help = "The number of additional terms to include in the fit (A,B,C,D). Must be in the range 0 to 4.")
    parser.add_argument('--no_plots', dest='no_plots', type=bool, required=False, default=False,
                        help = "Don't show plot of the results.")
    parser.add_argument('--outliers', dest='outliers', type=float, required=False, default=2.0,
                        help = "Sigma threshold for removing outliers in Wx, Wy.")

    args = parser.parse_args()    

    calibrate(args.hdf5, args.zoffsets, args.fit_order, args.outliers, args.no_plots)
   
