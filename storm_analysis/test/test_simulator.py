#!/usr/bin/env python
import numpy
import sys
import tifffile

import storm_analysis

import storm_analysis.spliner.spline_to_psf as splineToPSF

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate


def test_psf_spline2D_1():
    """
    Test that spline PSF agrees with spliner (for 0.0 offset).
    """
    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_name = storm_analysis.getData("test/data/test_spliner_psf_2d.spline")
        
    psf_sp_2d = psf.Spline2D(spline_name)
    sp_2d = splineToPSF.SplineToPSF2D(spline_name)

    psf_im = psf_sp_2d.getPSFPy(0.0, 0.0, 0.0)
    sp_im = sp_2d.getPSF(0.0, normalize = False)
    
    assert numpy.allclose(psf_im, sp_im)

    
def test_psf_spline2D_2():
    """
    Test that spline PSF C and Python versions agree.
    """
    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_name = storm_analysis.getData("test/data/test_spliner_psf_2d.spline")
        
    psf_sp_2D = psf.Spline2D(spline_name)
    dx = numpy.random.uniform(size = 5)
    dy = numpy.random.uniform(size = 5)

    for i in range(dx.size):
        psf_im_py = psf_sp_2D.getPSFPy(0.0, dy[i], dx[i])
        psf_im_c = psf_sp_2D.getPSF(0.0, dy[i], dx[i])

        assert numpy.allclose(psf_im_py, psf_im_c)


def test_psf_spline3D_1():
    """
    Test that spline PSF agrees with spliner (for 0.0 offset).
    """
    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_name = storm_analysis.getData("test/data/test_spliner_psf.spline")
        
    psf_sp_3d = psf.Spline3D(spline_name)
    sp_3d = splineToPSF.SplineToPSF3D(spline_name)

    psf_im = psf_sp_3d.getPSFPy(0.1, 0.0, 0.0)
    sp_im = sp_3d.getPSF(0.1, normalize = False)

    assert numpy.allclose(psf_im, sp_im)

        
def test_psf_spline3D_2():
    """
    Test that spline PSF C and Python versions agree.
    """
    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    spline_name = storm_analysis.getData("test/data/test_spliner_psf.spline")
        
    psf_sp_3D = psf.Spline3D(spline_name)
    dx = numpy.random.uniform(size = 5)
    dy = numpy.random.uniform(size = 5)
    zv = numpy.random.uniform(low = -0.5, high = 0.5, size = 5)

    for i in range(dx.size):
        psf_im_py = psf_sp_3D.getPSFPy(zv[i], dy[i], dx[i])
        psf_im_c = psf_sp_3D.getPSF(zv[i], dy[i], dx[i])

        assert numpy.allclose(psf_im_py, psf_im_c)
    
    
def test_simulate_1():
    """
    No photo-physics, simple PSF, ideal camera.
    """
    dax_name = storm_analysis.getPathOutputTest("test_sim1.dax")
    bin_name = storm_analysis.getData("test/data/test_sim.hdf5")

    sim = simulate.Simulate(background_factory = lambda settings, xs, ys, i3data : background.UniformBackground(settings, xs, ys, i3data),
                            camera_factory = lambda settings, xs, ys, i3data : camera.Ideal(settings, xs, ys, i3data, 100.0),
                            photophysics_factory = lambda settings, xs, ys, i3data : photophysics.AlwaysOn(settings, xs, ys, i3data, 1000.0),
                            psf_factory = lambda settings, xs, ys, i3data : psf.GaussianPSF(settings, xs, ys, i3data, 160.0),
                            x_size = 100, y_size = 75)

    sim.simulate(dax_name, bin_name, 5)


def test_simulate_2():
    """
    (Simple) STORM photo-physics, pure astigmatism PSF, EMCCD camera.
    """
    dax_name = storm_analysis.getPathOutputTest("test_sim2.dax")
    bin_name = storm_analysis.getData("test/data/test_sim.hdf5")

    sim = simulate.Simulate(background_factory = lambda settings, xs, ys, i3data : background.UniformBackground(settings, xs, ys, i3data, photons = 20),
                            camera_factory = lambda settings, xs, ys, i3data : camera.EMCCD(settings, xs, ys, i3data, 100.0, emccd_gain = 5.0, preamp_gain = 1.0, read_noise = 5),
                            photophysics_factory = lambda settings, xs, ys, i3data : photophysics.SimpleSTORM(settings, xs, ys, i3data, 4000.0, off_time = 10.0),
                            psf_factory = lambda settings, xs, ys, i3data : psf.PupilFunction(settings, xs, ys, i3data, 160.0, [[1.3, 2, 2]]),
                            x_size = 100, y_size = 75)
                   
    sim.simulate(dax_name, bin_name, 5)


def test_simulate_3():
    """
    No photo-physics, spline PSF, sCMOS camera.
    """
    
    # Only test for Python3 due to pickle incompatibility issues.
    if (sys.version_info < (3, 0)):
        return
    
    dax_name = storm_analysis.getPathOutputTest("test_sim3.dax")
    bin_name = storm_analysis.getData("test/data/test_sim.hdf5")
    cal_name = storm_analysis.getData("test/data/calib.npy")
    spline_name = storm_analysis.getData("test/data/test_spliner_psf.spline")

    sim = simulate.Simulate(background_factory = lambda settings, xs, ys, i3data : background.UniformBackground(settings, xs, ys, i3data, photons = 20),
                            camera_factory = lambda settings, xs, ys, i3data : camera.SCMOS(settings, xs, ys, i3data, cal_name),
                            photophysics_factory = lambda settings, xs, ys, i3data : photophysics.AlwaysOn(settings, xs, ys, i3data, 2000.0),
                            psf_factory = lambda settings, xs, ys, i3data : psf.Spline(settings, xs, ys, i3data, 160.0, spline_name))
                   
    sim.simulate(dax_name, bin_name, 5)

    
if (__name__ == "__main__"):
    test_psf_spline2D_1()
    test_psf_spline2D_2()
    test_psf_spline3D_1()
    test_psf_spline3D_2()
    test_simulate_1()
    test_simulate_2()
    test_simulate_3()
