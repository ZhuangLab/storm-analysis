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

import numpy

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_library.writeinsight3 as i3w


#
# A track stores most of its fields as the sum over the localizations in the
# track, see sa_utilities/tracker.py. Only x, y and z are averaged for you.
# These are the fields that have to be divided by the track length before they
# mean anything, i.e. the ones that describe a single observation of the
# molecule rather than a total over the track.
#
# 'height' and 'sum' are deliberately not in this list, they are reported as
# totals over the track. This is what the Insight3 era averager did. Its
# average_flag table in sa_utilities/avemlist.c, removed in b7470d4a, marked
# HEIGHT and AREA as TOTAL, and only WIDTH, ASPECT and BACKG as AVERAGE.
# Totaling the height is the older of the two, see 85a8fe90.
#
TRACK_MEAN_FIELDS = ["background", "error", "xsigma", "ysigma"]


def normalizeTrackFields(tracks):
    """
    Return a copy of 'tracks' with the per observation fields converted from
    sums to means. Without this a molecule that stayed on for ten frames is
    exported ten times too wide, sitting on ten times the background.
    """
    normalized = dict(tracks)
    if not "track_length" in tracks:
        return normalized

    length = tracks["track_length"].astype(numpy.float64)
    for field in TRACK_MEAN_FIELDS:
        if field in normalized:
            normalized[field] = normalized[field]/length

    return normalized


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
                i3.addMultiFitMolecules(normalizeTrackFields(tracks), 1, nm_per_pixel)

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
    
