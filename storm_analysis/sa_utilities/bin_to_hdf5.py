#!/usr/bin/env python
"""
Convert an Insight3 format binary file to a HDF5 storm-analysis format file.

This expects a mlist.bin file and not the tracked alist.bin file. Also the
Insight3 file must include metadata.

Hazen 12/17
"""
import os
import sys

from xml.etree import ElementTree

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.sa_h5py as saH5Py


def BinToHDF5(bin_name, hdf5_name):

    if os.path.exists(hdf5_name):
        os.remove(hdf5_name)

    with saH5Py.SAH5Py(hdf5_name, is_existing = False) as h5:

        if True:
            metadata = readinsight3.loadI3Metadata(bin_name)
            assert (metadata is not None), "No metadata found for " + bin_name
        
            movie_xml = metadata.find("movie")
            movie_x = int(movie_xml.find("movie_x").text)
            movie_y = int(movie_xml.find("movie_y").text)
            movie_l = int(movie_xml.find("movie_l").text)
            hash_value = movie_xml.find("hash_value")
            
            params_xml = metadata.find("settings")
            nm_per_pixel = float(params_xml.find("pixel_size").text)
            
        else:
            # Fill this in to work around missing metadata.
            movie_x = 256
            movie_y = 256
            movie_l = 10
            hash_value = ""

            params_xml = None
            nm_per_pixel = 100.0

        # Set movie properties.
        h5.setMovieInformation(movie_x, movie_y, movie_l, hash_value)

        # Set pixel size.
        h5.setPixelSize(nm_per_pixel)

        # Set metadata.
        if params_xml is not None:
            if (sys.version_info > (3, 0)):
                h5.addMetadata(ElementTree.tostring(params_xml, 'unicode'))
            else:
                h5.addMetadata(ElementTree.tostring(params_xml, 'ISO-8859-1'))

        # Convert data.
        i3 = readinsight3.I3Reader(bin_name)
        for i in range(i3.getNumberFrames()):
            i3data = i3.getMoleculesInFrame(i+1)
            h5_locs = i3dtype.convertToSAHDF5(i3data, i+1, nm_per_pixel)
            h5.addLocalizations(h5_locs, i)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Insight3 to HDF5 converter.')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The Insight3 file to convert.")
    parser.add_argument('--hdf5', dest='hdf5', type=str, required=True,
                        help = "The name for the HDF5 format file.")

    args = parser.parse_args()
    
    BinToHDF5(args.mlist, args.hdf5)
    
