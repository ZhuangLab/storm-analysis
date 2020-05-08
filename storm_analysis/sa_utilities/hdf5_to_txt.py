#!/usr/bin/env python
"""
Convert a HDF5 storm-analysis format file to csv text file.

If the HDF5 file has been tracked the text file will be created from
the tracks.

If the HDF5 file has not been tracked then the text file will be
created from the localizations.

Hazen 1/18
"""
import numpy
import os

from xml.etree import ElementTree

import storm_analysis.sa_library.sa_h5py as saH5Py


def hdf5ToTxt(hdf5_name, txt_name):
    with saH5Py.SAH5Reader(hdf5_name) as h5:
        nm_per_pixel = h5.getPixelSize()
        [movie_x, movie_y, movie_l, hash_value] = h5.getMovieInformation()

        with open(txt_name, "w") as fp:
            has_header = False
            fields = None
            
            # Convert tracks.
            if h5.hasTracks():
                index = 0
                print("Converting tracks.")
                for tracks in h5.tracksIterator():
                    
                    if not has_header:
                        fields = sorted(tracks.keys())
                        fp.write(",".join(["index"] + fields) + "\n")
                        has_header = True

                    for i in range(tracks["x"].size):
                        text = [str(index)]
                        for field in fields:
                            if(tracks[field].dtype == numpy.int32):
                                text.append(str(tracks[field][i]))
                            else:
                                text.append("{0:.3f}".format(tracks[field][i]))
                        fp.write(",".join(text) + "\n")
                        index += 1
                        
            # Convert localizations.
            else:
                index = 0
                print("Converting localizations.")
                for fnum, locs in h5.localizationsIterator(drift_corrected = False):
                                        
                    if not has_header:
                        fields = sorted(locs.keys())
                        fp.write(",".join(["index", "frame"] + fields) + "\n")
                        has_header = True

                    for i in range(locs["x"].size):
                        text = [str(index),str(fnum)]
                        for field in fields:
                            if(locs[field].dtype == numpy.int32):
                                text.append(str(locs[field][i]))
                            else:
                                text.append("{0:.3f}".format(locs[field][i]))
                        fp.write(",".join(text) + "\n")
                        index += 1


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'HDF5 to text converter.')

    parser.add_argument('--hdf5', dest='hdf5', type=str, required=True,
                        help = "The hdf5 file to convert.")
    parser.add_argument('--txt', dest='txt', type=str, required=True,
                        help = "The name for the text file.")

    args = parser.parse_args()
    
    hdf5ToTxt(args.hdf5, args.txt)
    
