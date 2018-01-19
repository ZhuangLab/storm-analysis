#!/usr/bin/env python
"""
Test 3D-DAOSTORM astigmatism z calibration functionality.
"""
import copy
import numpy

import storm_analysis

import storm_analysis.daostorm_3d.z_calibration as zCalibration
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_utilities.fitz_c as fitzC


def test_zcal_1():
    """
    Test simple fitting.
    """
    zv = numpy.arange(-0.6, 0.601, 0.01)
    z_params = [3.0, 0.3, 0.5]
    ww = zCalibration.zcalib0(z_params, zv)

    # No additional parameters.
    zf = zCalibration.doFit(ww, zv, copy.deepcopy(z_params), 0)
    assert(numpy.allclose(numpy.array(z_params), numpy.array(zf), atol = 0.01, rtol = 0.01))

    zf = zCalibration.doFit(ww, zv, copy.deepcopy(z_params), 2)
    z_params.append(0.0)
    z_params.append(0.0)
    assert(numpy.allclose(numpy.array(z_params), numpy.array(zf), atol = 0.01, rtol = 0.01))

    
def test_zcal_2():
    """
    Test X/Y swap handling.
    """
    zv = numpy.arange(-0.6, 0.601, 0.01)
    z_params = [3.0, 0.3, 0.5]
    ww = zCalibration.zcalib0(z_params, zv)

    # No additional parameters.
    inv_params = copy.deepcopy(z_params)
    inv_params[1] = -inv_params[1]
    zf = zCalibration.doFit(ww, zv, inv_params, 0)
    assert(numpy.allclose(numpy.array(z_params), numpy.array(zf), atol = 0.01, rtol = 0.01))


def test_zcal_3():
    """
    Test fitting both wx and wy.
    """
    zv = numpy.arange(-0.6, 0.601, 0.01)
    zx_params = [3.0, 0.3, 0.5]
    wx = zCalibration.zcalib0(zx_params, zv)
    zy_params = [3.0, -0.3, 0.5]
    wy = zCalibration.zcalib0(zy_params, zv)

    [wxp, wyp] = zCalibration.fitDefocusingCurves(wx, wy, zv)
    assert(numpy.allclose(numpy.array(zx_params), numpy.array(wxp), atol = 0.01, rtol = 0.01))
    assert(numpy.allclose(numpy.array(zy_params), numpy.array(wyp), atol = 0.01, rtol = 0.01))


def test_zcal_4():
    """
    Test measuring loading wx, wy from a file.
    """
    pixel_size = 100.0
    zv = numpy.arange(-0.6, 0.601, 0.01)
    zx_params = [3.0, 0.3, 0.5]
    wx = zCalibration.zcalib0(zx_params, zv)
    zy_params = [3.0, -0.3, 0.5]
    wy = zCalibration.zcalib0(zy_params, zv)

    # Write HDF5 file.
    h5_name = storm_analysis.getPathOutputTest("test_3ddao_zcal.hdf5")
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(100,100,wx.size,"")
        h5.setPixelSize(pixel_size)
        for i in range(wx.size):
            locs = {"x" : numpy.array([1.0]),
                    "y" : numpy.array([1.0]),
                    "xsigma" : numpy.array([0.5 * wx[i]]),
                    "ysigma" : numpy.array([0.5 * wy[i]])}
            h5.addLocalizations(locs, i)

    # Write offset file.
    off_name = storm_analysis.getPathOutputTest("test_3ddao_zcal.off")
    z_data = numpy.ones((wx.size, 2))
    z_data[:,1] = zv
    numpy.savetxt(off_name, z_data)

    # Load data.
    [rt_wx, rt_wy, rt_z, px_size] = zCalibration.loadWxWyZData(h5_name, off_name)

    assert(numpy.allclose(rt_wx, wx))
    assert(numpy.allclose(rt_wy, wy))
    assert(numpy.allclose(rt_z, zv))
    assert(px_size == pixel_size)

    # Write offset file with some 'bad' frames.
    off_name = storm_analysis.getPathOutputTest("test_3ddao_zcal.off")
    z_data = numpy.ones((wx.size, 2))
    z_data[:,1] = zv
    z_data[:10,0] = 0
    numpy.savetxt(off_name, z_data)

    # Load data.
    [rt_wx, rt_wy, rt_z, px_size] = zCalibration.loadWxWyZData(h5_name, off_name)
    assert (rt_z.size == (zv.size - 10))


def test_zcal_5():
    """
    Test that the z_calibration and fitz_c agree on the parameters.
    """
    pixel_size = 100.0
    zv = numpy.arange(-0.6, 0.601, 0.01)
    wx_params = [3.0, 0.3, 0.5]
    wx = zCalibration.zcalib0(wx_params, zv)
    wy_params = [3.0, -0.3, 0.5]
    wy = zCalibration.zcalib0(wy_params, zv)

    # Write HDF5 file.
    h5_name = storm_analysis.getPathOutputTest("test_3ddao_zcal.hdf5")
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(100, 100, wx.size, "")
        h5.setPixelSize(pixel_size)
        for i in range(wx.size):
            locs = {"x" : numpy.array([1.0]),
                    "y" : numpy.array([1.0]),
                    "xsigma" : numpy.array([0.5 * wx[i]]),
                    "ysigma" : numpy.array([0.5 * wy[i]])}
            h5.addLocalizations(locs, i)

    # Format parameters for fitzC.
    #
    wx_params = zCalibration.convertUnits(wx_params, pixel_size)
    wy_params = zCalibration.convertUnits(wy_params, pixel_size)

    # Fit z values.
    #
    fitzC.fitz(h5_name, 1.0, wx_params, wy_params, -0.7, 0.7)

    # Check results.
    fz = saH5Py.loadLocalizations(h5_name, fields = ["z"])["z"]

    assert(numpy.allclose(zv, fz, atol = 0.01, rtol = 0.01))


def test_zcal_6():
    """
    Test XML printing at least doesn't crash.
    """
    wx_params = [3.0, 0.3, 0.5]
    wy_params = [3.0, -0.3, 0.5]
    zCalibration.prettyPrint(wx_params, wy_params, 100.0)
    

if (__name__ == "__main__"):
    test_zcal_1()
    test_zcal_2()
    test_zcal_3()
    test_zcal_4()
    test_zcal_5()
    test_zcal_6()
    

    
