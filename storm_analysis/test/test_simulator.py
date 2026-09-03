#!/usr/bin/env python
import numpy
import sys
import tifffile

import storm_analysis

import storm_analysis.spliner.spline_to_psf as splineToPSF

import storm_analysis.simulator.astigmaticPSF as astigmaticPSF
import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.dhPSF as dhPSF
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

    
def test_astigmatic_psf():
    """
    astigmaticPSF takes z in microns, like everything else in the simulator.

    It used to multiply z by 0.001 before calling calcSxSy(), which says the
    incoming z is nanometers. Fed the simulator's own convention that put
    every z within a nanometer of focus, so it returned a round,
    z-independent PSF, sx = sy = 1.068 px at every z.
    """
    z = numpy.array([-0.4, -0.2, 0.0, 0.2, 0.4])
    x = numpy.ones(z.size)
    y = numpy.ones(z.size)
    h = numpy.ones(z.size)*100.0

    objects = astigmaticPSF.PSF(x, y, z, h)
    sx = objects[:,3]
    sy = objects[:,4]

    # Round at focus.
    assert(abs(sx[2] - sy[2]) < 1.0e-9)

    # Elongated in x below focus, in y above it, and it is a real effect and
    # not rounding. wx_params puts the x waist at +150nm and the y waist at
    # -150nm, with a 400nm defocusing scale.
    assert(sx[0] > 1.4*sy[0])
    assert(sy[4] > 1.4*sx[4])

    # Symmetric about focus.
    assert(abs(sx[0] - sy[4]) < 1.0e-9)
    assert(abs(sy[0] - sx[4]) < 1.0e-9)

    # And the width actually varies with z, which is what a Z fitter needs.
    assert(sx[0] > 1.5*sx[2])

    # PSFIntegral uses the same widths.
    integral = astigmaticPSF.PSFIntegral(z, h)
    assert(numpy.allclose(integral, 2.0*numpy.pi*h*sx*sy))


def test_dh_psf():
    """
    dhPSF takes z in microns too.

    z_min and z_max were -500 and 500, which are nanometers, so feeding it
    the simulator's own z froze the helix angle at 45 degrees across the
    whole range. A double helix PSF that does not rotate carries no z
    information at all, the same degeneracy astigmaticPSF had.
    """
    [xc, yc] = [10.0, 20.0]
    z = numpy.array([-0.5, 0.0, 0.5])

    objects = dhPSF.PSF(numpy.ones(z.size)*xc,
                        numpy.ones(z.size)*yc,
                        z,
                        numpy.ones(z.size)*100.0)

    # Two lobes per emitter, ordered [emitter0 +, emitter0 -, emitter1 +, ..].
    assert(objects.shape[0] == 2*z.size)

    lobe = objects[::2,:]
    dx = lobe[:,0] - xc
    dy = lobe[:,1] - yc

    r = 2.0*dhPSF.sigma

    # At z_min the lobes lie along x, at z_max along y, and at focus the
    # angle is halfway between.
    assert(numpy.allclose([dx[0], dy[0]], [r, 0.0], atol = 1.0e-9))
    assert(numpy.allclose([dx[2], dy[2]], [0.0, r], atol = 1.0e-9))
    assert(abs(dx[1] - dy[1]) < 1.0e-9)

    # The angle has to actually sweep, which is the whole point.
    angles = numpy.degrees(numpy.arctan2(dy, dx))
    assert(numpy.allclose(angles, [0.0, 45.0, 90.0], atol = 1.0e-6))

    # The second lobe of each pair is opposite the first.
    other = objects[1::2,:]
    assert(numpy.allclose(other[:,0] - xc, -dx))
    assert(numpy.allclose(other[:,1] - yc, -dy))


if (__name__ == "__main__"):
    test_psf_spline2D_1()
    test_psf_spline2D_2()
    test_psf_spline3D_1()
    test_psf_spline3D_2()
    test_simulate_1()
    test_simulate_2()
    test_simulate_3()
    test_astigmatic_psf()
    test_dh_psf()
