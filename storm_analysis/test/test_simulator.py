#!/usr/bin/env python

import storm_analysis

import storm_analysis.simulator.background as background
import storm_analysis.simulator.camera as camera
import storm_analysis.simulator.photophysics as photophysics
import storm_analysis.simulator.psf as psf
import storm_analysis.simulator.simulate as simulate

def test_simulate_1():
    """
    No photo-physics, simple PSF, ideal camera.
    """
    dax_name = storm_analysis.getPathOutputTest("test_sim1.dax")
    bin_name = storm_analysis.getData("test/data/test_sim_olist.bin")

    sim = simulate.Simulate(lambda settings, xs, ys, i3data : background.UniformBackground(settings, xs, ys, i3data),
                            lambda settings, xs, ys, i3data : camera.Ideal(settings, xs, ys, i3data, 100.0),
                            lambda settings, xs, ys, i3data : photophysics.AlwaysOn(settings, xs, ys, i3data, 1000.0),
                            lambda settings, xs, ys, i3data : psf.GaussianPSF(settings, xs, ys, i3data, 160.0),
                            x_size = 100, y_size = 75)
                   
    sim.simulate(dax_name, bin_name, 5)

def test_simulate_2():
    """
    (Simple) STORM photo-physics, pure astigmatism PSF, EMCCD camera.
    """
    dax_name = storm_analysis.getPathOutputTest("test_sim2.dax")
    bin_name = storm_analysis.getData("test/data/test_sim_olist.bin")

    sim = simulate.Simulate(lambda settings, xs, ys, i3data : background.UniformBackground(settings, xs, ys, i3data, photons = 20),
                            lambda settings, xs, ys, i3data : camera.EMCCD(settings, xs, ys, i3data, 100.0, emccd_gain = 5.0, preamp_gain = 1.0, read_noise = 5),
                            lambda settings, xs, ys, i3data : photophysics.SimpleSTORM(settings, xs, ys, i3data, 4000.0, off_time = 10.0),
                            lambda settings, xs, ys, i3data : psf.PupilFunction(settings, xs, ys, i3data, 160.0, [[1.3, 2, 2]]),
                            x_size = 100, y_size = 75)
                   
    sim.simulate(dax_name, bin_name, 5)

def test_simulate_3():
    """
    No photo-physics, spline PSF, sCMOS camera.
    """
    dax_name = storm_analysis.getPathOutputTest("test_sim3.dax")
    bin_name = storm_analysis.getData("test/data/test_sim_olist.bin")
    cal_name = storm_analysis.getData("test/data/calib.npy")
    spline_name = storm_analysis.getData("test/data/test_spliner_psf.spline")

    sim = simulate.Simulate(lambda settings, xs, ys, i3data : background.UniformBackground(settings, xs, ys, i3data, photons = 20),
                            lambda settings, xs, ys, i3data : camera.SCMOS(settings, xs, ys, i3data, 100.0, cal_name),
                            lambda settings, xs, ys, i3data : photophysics.AlwaysOn(settings, xs, ys, i3data, 2000.0),
                            lambda settings, xs, ys, i3data : psf.Spline(settings, xs, ys, i3data, 160.0, spline_name))
                   
    sim.simulate(dax_name, bin_name, 5)

    


    
if (__name__ == "__main__"):
    test_simulate_1()
    test_simulate_2()
    test_simulate_3()

    
