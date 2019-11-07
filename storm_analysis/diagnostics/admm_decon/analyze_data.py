#!/usr/bin/env python
"""
Analyze test data using ADMM decon. This just runs ADMM decon
without any fitting as ADMM isn't currently part of any analysis
pipeline. It is also limited to 2D.

Hazen 10/19
"""
import glob
import numpy
import os
import tifffile
import time

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.parameters as params

import storm_analysis.spliner.spline_to_psf as splineToPSF

import storm_analysis.rolling_ball_bgr.rolling_ball as rollingBall
    
import storm_analysis.admm.admm_decon as admmDecon


def analyze(movie_name, mlist_name, settings_name):
    """
    Simple wrapper for testing ADMM decon.
    """
    # Load parameters
    parameters = params.ParametersSplinerFISTA().initFromFile(settings_name)

    # Load movie.
    frame_reader = analysisIO.FrameReaderStd(movie_file = movie_name,
                                             parameters = parameters)
    movie_reader = analysisIO.MovieReader(frame_reader = frame_reader,
                                          parameters = parameters)
    image = frame_reader.loadAFrame(0)

    # Load spline.
    psf_object = splineToPSF.loadSpline(parameters.getAttr("spline"))

    # Create Rolling ball background removal object.
    rb = rollingBall.RollingBall(parameters.getAttr("rb_radius"),
                                 parameters.getAttr("rb_sigma"))
    
    # Create ADMM decon object.
    adecon = admmDecon.ADMMDecon(image.shape,
                                 psf_object,
                                 parameters.getAttr("fista_number_z"),
                                 parameters.getAttr("fista_timestep"))

    # Create HDF5 data writer.
    data_writer = analysisIO.DataWriterHDF5(data_file = mlist_name,
                                            parameters = parameters,
                                            sa_type = "ADMM")

    # Analyze movie.
    movie_reader.setup(-1)
    with tifffile.TiffWriter(os.path.join(os.path.dirname(movie_name), "admm.tif")) as tf:
        while movie_reader.nextFrame():
            image = movie_reader.getFrame()

            # Estimate background.
            background = rb.estimateBG(image)

            adecon.newImage(image)
            adecon.newBackground(background)

            # Do deconvolution.
            adecon.decon(parameters.getAttr("fista_iterations"),
                         parameters.getAttr("fista_lambda"),
                         verbose = False)

            # For debugging.
            fx = adecon.getXVector()
            for i in range(fx.shape[2]):
                tf.save(fx[:,:,i].astype(numpy.float32))
                    
            # Save peaks.
            peaks = adecon.getPeaks(parameters.getAttr("fista_threshold"), 5)
            data_writer.addPeaks(peaks, movie_reader)

            # Feedback.
            print("Frame:", movie_reader.getCurrentFrameNumber(), data_writer.getNumberAdded(), data_writer.getTotalPeaks())

    # Close IO files.
    movie_reader.close()
    data_writer.close(True)
    
    # Clean up.
    adecon.cleanup()
    rb.cleanup()
    

def analyzeData():
    dirs = sorted(glob.glob("test*"))

    total_time = 0.0
    for a_dir in dirs:
        print()
        print("Analyzing:", a_dir)
        print()
    
        mlist = a_dir + "/test.hdf5"

        # Remove stale results, if any.
        if os.path.exists(mlist):
            os.remove(mlist)

        # Run analysis.
        start_time = time.time()
        analyze(a_dir + "/test.tif", mlist, "adecon.xml")
        stop_time = time.time()

        # Save timing results.
        total_time += stop_time - start_time
        print("Analysis completed in {0:.2f} seconds".format(stop_time - start_time))
        
        with open(a_dir + "/timing.txt", "w") as fp:
            fp.write(str(stop_time - start_time) + "\n")

    print()
    print("{0:d} directories analyzed in {1:.2f} seconds.".format(len(dirs), total_time))


if (__name__ == "__main__"):
    analyzeData()
    
