#!/usr/bin/env python
"""
Tests of multi_plane.fitting_mp
"""
import glob
import numpy
import os
import pickle
import tifffile

import storm_analysis
import storm_analysis.sa_library.fitting as fitting
import storm_analysis.sa_library.parameters as parameters

import storm_analysis.multi_plane.fitting_mp as fittingMp


im_size = (20, 18)
n_channels = 2


def configureTest():
    """
    Build the parameters and the mapping file that MPPeakFinderDao needs.
    """
    map_name = storm_analysis.getPathOutputTest("test_fitting_mp.map")

    mappings = {}
    for i in range(n_channels):
        mappings["0_" + str(i) + "_x"] = numpy.array([0.0, 1.0, 0.0])
        mappings["0_" + str(i) + "_y"] = numpy.array([0.0, 0.0, 1.0])
        mappings[str(i) + "_0_x"] = numpy.array([0.0, 1.0, 0.0])
        mappings[str(i) + "_0_y"] = numpy.array([0.0, 0.0, 1.0])
    with open(map_name, "wb") as fp:
        pickle.dump(mappings, fp)

    params = parameters.ParametersMultiplaneDao()
    params.setAttr("background_sigma", "float", 8.0)
    params.setAttr("find_max_radius", "int", 2)
    params.setAttr("foreground_sigma", "float", 1.5)
    params.setAttr("iterations", "int", 1)
    params.setAttr("mapping", "filename", map_name)
    params.setAttr("no_fitting", "int", 0)
    params.setAttr("pixel_size", "float", 100.0)
    params.setAttr("roi_size", "int", 10)
    params.setAttr("sigma", "float", 1.5)
    params.setAttr("threshold", "float", 6.0)

    return params


def test_fitting_mp_check_mode():
    """
    MPPeakFinderDao.setVariances() saves a picture of each channel's PSF
    when check_mode is on.

    That block was copied from the arbitrary PSF version, where the PSF is
    called 'psf' and the channel index is 'j'. Here the PSF is 'psf_norm'
    and the index is 'i', so it raised NameError, and the filename kept
    "{1:d}" from the two argument version, so it would have raised
    IndexError once the names were fixed.
    """
    params = configureTest()

    # check_mode writes into the working directory.
    cwd = os.getcwd()
    os.chdir(storm_analysis.getPathOutputTest(""))
    try:
        for name in glob.glob("psf_z0.0_c*.tif"):
            os.remove(name)

        finder = fittingMp.MPPeakFinderDao(parameters = params,
                                           n_channels = n_channels)
        finder.check_mode = True

        variances = [numpy.ones(im_size) for i in range(n_channels)]
        finder.setVariances(variances)

        # One image per channel, named by the channel index.
        written = sorted(glob.glob("psf_z0.0_c*.tif"))
        assert(written == ["psf_z0.0_c0.tif", "psf_z0.0_c1.tif"])

        # And it is the normalized PSF that was saved, computed here
        # independently of the finder. The variances are padded by the
        # margin before the filters are built.
        padded = (im_size[0] + 2*finder.margin, im_size[1] + 2*finder.margin)
        expected = fitting.gaussianPSF(padded, 1.5).astype(numpy.float32)

        for name in written:
            psf = tifffile.imread(name)
            assert(psf.shape == padded)
            assert(numpy.allclose(psf, expected))
    finally:
        os.chdir(cwd)


if (__name__ == "__main__"):
    test_fitting_mp_check_mode()
