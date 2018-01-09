#!/usr/bin/env python
"""
This is mostly for debugging. It takes the original HDF5
localization file and makes one for each channel.

Hazen 01/18
"""
import numpy
import os

import storm_analysis.sa_library.sa_h5py as saH5Py


def separateChannels(h5_name):

    h5w = []
    with saH5Py.SAH5Py(h5_name) as h5:

        # Create a writer for each channel.
        assert(h5.getNChannels() > 1), "Data only has a single channel."
        for i in range(h5.getNChannels()):
            temp = h5_name[:-5] + "_c" + str(i) + ".hdf5"
            h5w_temp = saH5Py.SAH5Py(temp, is_existing = False, overwrite = True)
            h5w_temp.setPixelSize(h5.getPixelSize())
            h5w_temp.setMovieInformation(*h5.getMovieInformation())
            h5w.append(h5w_temp)

        # Split out data for each channel.
        for fnum, locs in h5.localizationsIterator(drift_corrected = False):
            split_locs = h5.splitByChannel(locs)
            
            for i in range(len(h5w)):
                h5w[i].addLocalizations(split_locs[i], fnum)
                
    # Close writers.
    for elt in h5w:
        elt.close()

        
if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Separate out the channels of a HDF5 file.')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the storm-analysis HDF5 file.")

    args = parser.parse_args()
    
    separateChannels(args.mlist)
