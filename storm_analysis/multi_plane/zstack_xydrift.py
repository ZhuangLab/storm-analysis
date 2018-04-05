#!/usr/bin/env python
"""
Estimate the xy drift that occured during a (PSF measurement)
z stack. This just determines the average position difference
between the beads in the first frame and the beads in the last
frame. The beads need to be more or less in focus in the first
and last frame. Also we assume that the total drift is less
than a single pixel.

Hazen 09/17
"""
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_library.ia_utilities_c as utilC

def xyDrift(locs_filename):

    # Load localizations.
    #
    with saH5Py.SAH5Py(locs_filename) as h5:
        n_frames = h5.getMovieLength()
        f0_locs = h5.getLocalizationsInFrame(0, fields = ["x", "y"])
        fn_locs = h5.getLocalizationsInFrame(n_frames - 1, fields = ["x", "y"])

    assert (f0_locs["x"].size > 0), "No localizations in the first frame."
    assert (fn_locs["y"].size > 0), "No localizations in the last frame."

    #
    # Identify matching beads in the first and last frame and
    # compute the displacement.
    #
    p_index = utilC.peakToPeakDistAndIndex(f0_locs['x'], f0_locs['y'],
                                           fn_locs['x'], fn_locs['y'])[1]

    all_dx = []
    all_dy = []
    for i in range(f0_locs['x'].size):
        dx = fn_locs['x'][p_index[i]] - f0_locs['x'][i]
        dy = fn_locs['y'][p_index[i]] - f0_locs['y'][i]
        if ((dx*dx + dy*dy) < 1.0):
            all_dx.append(dx)
            all_dy.append(dy)

    #
    # Return average per frame.
    #
    return [numpy.mean(numpy.array(all_dx))/float(n_frames),
            numpy.mean(numpy.array(all_dy))/float(n_frames)]


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Estimate XY drift in during a PSF measurement bead stack.')

    parser.add_argument('--bin', dest='in_bin', type=str, required=True)

    args = parser.parse_args()

    [dx, dy] = xyDrift(args.in_bin)
    print("dx: {0:5f} dy: {1:5f}".format(dx, dy))

