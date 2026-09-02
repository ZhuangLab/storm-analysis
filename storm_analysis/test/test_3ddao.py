#!/usr/bin/env python
import numpy

import storm_analysis

import storm_analysis.daostorm_3d.find_peaks as findPeaks
import storm_analysis.sa_library.dao_fit_c as daoFitC
import storm_analysis.sa_library.parameters as params

import storm_analysis.test.verifications as veri


def test_3ddao_2d_fixed():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_fixed.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1998):
        raise Exception("3D-DAOSTORM 2D fixed did not find the expected number of localizations.")
    

def test_3ddao_2d_fixed_gt():
    """
    Start fitting from ground truth locations.
    """
    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed_gt.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_fixed_gt.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 200):
        raise Exception("3D-DAOSTORM 2D fixed ground truth did not find the expected number of localizations.")


def test_3ddao_2d_fixed_gt_text():
    """
    Start fitting from ground truth locations (text file version).
    """
    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed_gt_text.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_fixed_gt_text.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 200):
        raise Exception("3D-DAOSTORM 2D fixed ground truth did not find the expected number of localizations.")    
    
    
def test_3ddao_2d_fixed_low_snr():

    movie_name = storm_analysis.getData("test/data/test_low_snr.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed_low_snr.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_fixed_low_snr.hdf5")
    storm_analysis.removeFile(mlist)
    
    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 392):
        raise Exception("3D-DAOSTORM 2D fixed low snr did not find the expected number of localizations.")
    

def test_3ddao_2d_fixed_non_square():

    movie_name = storm_analysis.getData("test/data/test_300x200.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d_300x200.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1002):
        raise Exception("3D-DAOSTORM 2D fixed non square did not find the expected number of localizations.")
    
    
def test_3ddao_2d():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_2d.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_2d.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1970):
        raise Exception("3D-DAOSTORM 2D did not find the expected number of localizations.")


def test_3ddao_3d():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_3d.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_3d.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1956):
        raise Exception("3D-DAOSTORM 3D did not find the expected number of localizations.")    


def test_3ddao_Z():

    movie_name = storm_analysis.getData("test/data/test.dax")
    settings = storm_analysis.getData("test/data/test_3d_Z.xml")
    mlist = storm_analysis.getPathOutputTest("test_3d_Z.hdf5")
    storm_analysis.removeFile(mlist)

    from storm_analysis.daostorm_3d.mufit_analysis import analyze
    analyze(movie_name, mlist, settings)

    # Verify number of localizations found.
    num_locs = veri.verifyNumberLocalizations(mlist)
    if not veri.verifyIsCloseEnough(num_locs, 1955):
        raise Exception("3D-DAOSTORM Z did not find the expected number of localizations.")


def test_3ddao_Z_jacobian():
    """
    Check the 'Z' model's analytic derivative of the error with respect to z
    against a numerical derivative of the same quantity.

    The 'Z' model is the only one where z is a fitting parameter, so it is the
    only place where the width versus z calibration gets differentiated. An
    error there does not move the converged answer, since the fit stops where
    the residual is minimized either way, it just costs iterations. That makes
    it invisible to the localization count tests, hence this one.
    """
    # Deliberately different defocusing terms in x and y so that an error in
    # either one shows up.
    wx_params = numpy.array([3.34, -0.30, 0.40, 0.0, 0.0])
    wy_params = numpy.array([3.34, 0.30, 0.55, 0.0, 0.0])

    im_size = 40
    [x, y] = [20.0, 20.0]
    [height, background] = [2000.0, 20.0]

    def sigmaFromZ(w_params, z):
        zt = (z - w_params[1])/w_params[2]
        tmp = 1.0 + zt*zt + w_params[3]*zt*zt*zt + w_params[4]*zt*zt*zt*zt
        return 0.5*w_params[0]*numpy.sqrt(tmp)

    # Astigmatic Gaussian, using the same width versus z model as the fitter.
    [sx, sy] = [sigmaFromZ(wx_params, 0.05), sigmaFromZ(wy_params, 0.05)]
    [yi, xi] = numpy.mgrid[0:im_size,0:im_size]
    image = height*numpy.exp(-((xi-x)*(xi-x)/(2.0*sx*sx) + (yi-y)*(yi-y)/(2.0*sy*sy))) + background

    def errorAtZ(z):
        """
        Add a peak at 'z' and return both the fit error and the analytic
        derivative of the error with respect to z.
        """
        mfitter = daoFitC.MultiFitterZ(roi_size = 20,
                                       wx_params = wx_params,
                                       wy_params = wy_params,
                                       min_z = -0.6,
                                       max_z = 0.6,
                                       rqe = numpy.ones(image.shape),
                                       scmos_cal = numpy.zeros(image.shape))
        mfitter.initializeC(image)
        mfitter.newImage(image)
        mfitter.newBackground(numpy.ones(image.shape)*background)
        mfitter.newPeaks({"x" : numpy.array([x]),
                          "y" : numpy.array([y]),
                          "z" : numpy.array([z]),
                          "background" : numpy.array([background]),
                          "height" : numpy.array([height]),
                          "xsigma" : numpy.array([sigmaFromZ(wx_params, z)]),
                          "ysigma" : numpy.array([sigmaFromZ(wy_params, z)])},
                         "text")
        err = mfitter.getPeakProperty("error")[0]
        d_err = mfitter.getPeakProperty("jacobian")[0][3]
        mfitter.cleanup(verbose = False)
        return [err, d_err]

    # Evaluate away from the best fit z so that the derivative is not zero.
    z = 0.15
    dz = 1.0e-5

    analytic = errorAtZ(z)[1]
    numerical = (errorAtZ(z + dz)[0] - errorAtZ(z - dz)[0])/(2.0*dz)

    assert(abs(analytic/numerical - 1.0) < 1.0e-3)


def test_3ddao_scmos_cal():
    """
    Test that scmos calibration data is initialized to 0.0.
    """
    settings = storm_analysis.getData("test/data/test_3d_2d_fixed.xml")
    parameters = params.ParametersDAO().initFromFile(settings)

    # Create analysis object and reach deep into it..
    find_fit = findPeaks.initFindAndFit(parameters)
    fitter = find_fit.peak_fitter
    mfitter = fitter.mfitter

    # Initialize with an image.
    image = numpy.ones((100,100))
    fitter.newImage(image)

    # Verify that the image is still all ones.
    resp = mfitter.getResidual()
    assert(numpy.max(resp - 1.0) < 1.0e-6)

    # Cleanup.
    fitter.cleanUp()

    
if (__name__ == "__main__"):
    test_3ddao_2d_fixed()
    test_3ddao_2d_fixed_gt()
    test_3ddao_2d_fixed_gt_text()
    test_3ddao_2d_fixed_low_snr()
    test_3ddao_2d_fixed_non_square()
    test_3ddao_2d()
    test_3ddao_3d()
    test_3ddao_Z()
    test_3ddao_Z_jacobian()
    test_3ddao_scmos_cal()

