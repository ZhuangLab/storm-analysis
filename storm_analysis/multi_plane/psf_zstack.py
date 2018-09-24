#!/usr/bin/env python
"""
Given a movie and list of locations (the output of
multi_plane.psf_localizations), generate a list of
z-stacks.

The z stack results are in units of photo-electrons.

FIXME: Averaging should be done with weighting by pixel 
       variance?

FIXME: Drift correction, if specified, is not corrected for 
       the channel to channel mapping.

Hazen 02/18
"""
import numpy
import scipy
import scipy.ndimage
import tifffile

import storm_analysis.sa_library.analysis_io as analysisIO
import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.spliner.measure_psf_utils as measurePSFUtils


def psfZStack(movie_name, h5_filename, zstack_name, scmos_cal = None, aoi_size = 8, driftx = 0.0, drifty = 0.0):
    """
    movie_name - The movie file containing the z stack.
    h5_filename - The HDF5 file containing the localizations to use for the PSF measurement.
    zstack_name - The name of the file to save the zstack in.
    scmos_cal - The sCMOS calibration file.
    aoi_size - The AOI size in pixels.

    driftx, drifty are in units of pixels per frame, (bead x last frame - bead x first frame)/n_frames.
    """
    # Create appropriate reader.
    if scmos_cal is None:
        fr_reader = datareader.inferReader(movie_name)
    else:
        fr_reader = analysisIO.FrameReaderSCMOS(movie_file = movie_name,
                                                calibration_file = scmos_cal)
        
    [movie_x, movie_y, movie_len] = fr_reader.filmSize()
    
    # Load localizations.
    with saH5Py.SAH5Py(h5_filename) as h5:
        locs = h5.getLocalizations()
        x = locs["y"] + 1
        y = locs["x"] + 1

    # Measure Z stacks.
    z_stacks = []
    for i in range(x.size):
        z_stacks.append(numpy.zeros((2*aoi_size, 2*aoi_size, movie_len)))
        
    for i in range(movie_len):
        if((i%50)==0):
            print("Processing frame {0:0d}".format(i))

        # Load the frame. This also handles gain and offset correction.
        #
        frame = fr_reader.loadAFrame(i)

        # Subtract estimated background. This assumes that the image is
        # mostly background and that the background is uniform.
        #
        frame = frame - numpy.median(frame)
            
        for j in range(x.size):
            xf = x[j] + driftx * float(i)
            yf = y[j] + drifty * float(i)
            z_stacks[j][:,:,i] = measurePSFUtils.extractAOI(frame, aoi_size, xf, yf)

    # Save z_stacks.
    numpy.save(zstack_name + ".npy", z_stacks)

    # Save a (normalized) z_stack as tif for inspection purposes.
    z_stack = z_stacks[0]
    z_stack = z_stack/numpy.amax(z_stack)
    z_stack = z_stack.astype(numpy.float32)
    with tifffile.TiffWriter(zstack_name + ".tif") as tf:
        for i in range(movie_len):
            tf.save(z_stack[:,:,i])

            
if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Average AOIs together into a z-stack for PSF measurement.')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie, can be .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations psf file.")
    parser.add_argument('--zstack', dest='zstack', type=str, required=True,
                        help = "The name of the file to save the zstack (without an extension).")
    parser.add_argument('--scmos_cal', dest='scmos_cal', type=str, required=False, default = None,
                        help = "The name of the sCMOS calibration data file.")    
    parser.add_argument('--aoi_size', dest='aoi_size', type=int, required=False, default=8,
                        help = "The size of the area of interest around the bead in pixels. The default is 8.")
    parser.add_argument('--driftx', dest='driftx', type=float, required=False, default=0.0,
                        help = "Drift in x in pixels per frame. The default is 0.0.")
    parser.add_argument('--drifty', dest='drifty', type=float, required=False, default=0.0,
                        help = "Drift in y in pixels per frame. The default is 0.0.")

    args = parser.parse_args()
    
    psfZStack(args.movie,
              args.mlist,
              args.zstack,
              scmos_cal = args.scmos_cal,
              aoi_size = args.aoi_size,
              driftx = args.driftx,
              drifty = args.drifty)
