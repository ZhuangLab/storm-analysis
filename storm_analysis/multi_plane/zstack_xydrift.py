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

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.ia_utilities_c as utilC

def xyDrift(locs_filename):
    locs_i3 = readinsight3.I3Reader(locs_filename)

    #
    # Load the number of frames and the localizations in
    # the first and last frame.
    #
    n_frames = locs_i3.getNumberFrames()
    
    f0_locs = locs_i3.getMoleculesInFrame(1)
    fn_locs = locs_i3.getMoleculesInFrame(n_frames)

    assert (f0_locs.size > 0), "No localizations in the first frame."
    assert (fn_locs.size > 0), "No localizations in the last frame."

    #
    # Identify matching beads in the first and last frame and
    # compute the displacement.
    #
    p_index = utilC.peakToPeakIndex(f0_locs['xc'], f0_locs['yc'],
                                    fn_locs['xc'], fn_locs['yc'])

    all_dx = []
    all_dy = []
    for i in range(f0_locs.size):
        dx = fn_locs['xc'][p_index[i]] - f0_locs['xc'][i]
        dy = fn_locs['yc'][p_index[i]] - f0_locs['yc'][i]
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

