#!/usr/bin/env python
"""
Convert a HDF5 storm-analysis format file to an Insight3 format bin file.

If the HDF5 file has been tracked the Insight3 file will be created from
the tracks. Note however that in the converted Insight3 file all the 
localizations will be in frame 1.

If the HDF5 file has not been tracked then the Insight3 file will be
created from the localizations.

Hazen 12/17
"""
import os

from xml.etree import ElementTree

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_library.writeinsight3 as i3w


def hdf5ToBin(hdf5_name, bin_name):
    with saH5Py.SAH5Reader(hdf5_name) as h5:
        nm_per_pixel = h5.getPixelSize()
        [movie_x, movie_y, movie_l, hash_value] = h5.getMovieInformation()

        # Create Insight3 file for writing.
        i3 = i3w.I3Writer(bin_name)

        # Convert tracks.
        if h5.hasTracks():
            print("Converting tracks.")
            for tracks in h5.tracksIterator():
                i3.addMultiFitMolecules(tracks, 1, nm_per_pixel)

        # Convert localizations.
        else:
            print("Converting localizations.")
            for fnum, locs in h5.localizationsIterator(drift_corrected = False):
                i3.addMultiFitMolecules(locs, fnum + 1, nm_per_pixel)

        # Add metadata.
        etree = ElementTree.Element("xml")

        # Analysis parameters.
        h5_metadata = h5.getMetadata()
        etree.append(ElementTree.fromstring(h5_metadata))

        # Movie properties.
        movie_props = ElementTree.SubElement(etree, "movie")
        field = ElementTree.SubElement(movie_props, "hash_value")
        field.text = hash_value
        for elt in [["movie_x", movie_x],
                    ["movie_y", movie_y],
                    ["movie_l", movie_l]]:
            field = ElementTree.SubElement(movie_props, elt[0])
            field.text = str(elt[1])
            
        metadata = ElementTree.tostring(etree, 'ISO-8859-1')
        
        # Close i3 file with metadata.
        i3.closeWithMetadata(metadata)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'HDF5 to Insight3 converter.')

    parser.add_argument('--hdf5', dest='hdf5', type=str, required=True,
                        help = "The hdf5 file to convert.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name for the Insight3 format binary file.")

    args = parser.parse_args()
    
    hdf5ToBin(args.hdf5, args.mlist)
    
