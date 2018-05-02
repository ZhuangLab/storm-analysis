#!/usr/bin/env python
"""
Convert an HDF5 file to a beads file.

This is useful for example when you analyze a frame of a movie
with sCMOS or 3D-DAOSTORM and you want to use those localizations
for spline determination.

Hazen 05/18
"""
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py

def hdf5ToBeads(hdf5_file, beads_file, frame = 0):
    """
    hdf5_file - The HDF5 file containing the analysis results.
    beads_file - The beads.txt file to save the localizations in.
    frame - Which frame to get the localizations from.
    """
    with saH5Py.SAH5Py(hdf5_file) as h5:
        locs = h5.getLocalizationsInFrame(frame)

    assert(locs['x'].size > 0), "No localizations in frame " + str(frame)
    
    numpy.savetxt(beads_file,
                  numpy.transpose(numpy.vstack((locs['x'],
                                                locs['y'],
                                                locs['height'],
                                                locs['background']))))

    
if (__name__ == "__main__"):

    import argparse
    
    parser = argparse.ArgumentParser(description = 'Convert an HDF5 file to a beads.txt file for PSF measurement.')

    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the HDF5 binary file.")
    parser.add_argument('--beads', dest='beads', type=str, required=True,
                        help = "The name of the beads file.")
    parser.add_argument('--frame', dest='frame', type=int, required=False, default = 0,
                        help = "Which frame to use the localizations from. The default is 0.")

    args = parser.parse_args()

    hdf5ToBeads(args.hdf5,
                args.beads,
                frame = args.frame)
    
