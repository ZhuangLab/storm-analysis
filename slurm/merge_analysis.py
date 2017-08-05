#!/usr/bin/env python
"""
Merges the intermediate localization files into a single 
localization file.

Hazen 08/17
"""

import glob
import os
from xml.etree import ElementTree

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

def mergeAnalysis(dir_name, bin_name, ext):

    # Create Insight3 file writer.
    i3_out = writeinsight3.I3Writer(bin_name)

    # Find all the job*.xml files.
    job_xml_files = glob.glob(dir_name + "job*.xml")

    # Sort job files.
    job_xml_files = sorted(job_xml_files, key = lambda x: int(os.path.splitext(os.path.basename(x))[0].split("_")[1]))

    # Check for corresponding mlist.bin files.
    metadata = None
    last_frame = 0
    for i in range(len(job_xml_files)):

        mlist_name = dir_name + "p_" + str(i+1) + "_mlist" + ext

        if os.path.exists(mlist_name):
            if readinsight3.checkStatus(mlist_name):

                # Load metadata from the first file.
                if (i == 0):
                    metadata = readinsight3.loadI3Metadata(mlist_name)

                # Read localizations.
                i3_data = readinsight3.loadI3File(mlist_name, verbose = False)

                # Print frame range covered.
                last_frame = i3_data["fr"][-1]
                print(i3_data["fr"][0], last_frame, mlist_name)

                # Add localizations to the output file.
                i3_out.addMolecules(i3_data)

                continue

        print("Merge failed because", job_xml_files[i], "is incomplete.")
        i3_out.close()
        os.remove(bin_name)
        assert(False)

    if metadata is None:
        print("No metadata found.")
        i3_out.close()
    else:

        # Fix movie length node based on the last frame of the last molecule.
        metadata.find("movie").find("movie_l").text = str(last_frame)

        # Also need to fix analysis end points. We are assuming that the
        # entire movie was analyzed.
        metadata.find("settings").find("start_frame").text = "-1"
        metadata.find("settings").find("max_frame").text = "-1"

        i3_out.closeWithMetadata(ElementTree.tostring(metadata, 'ISO-8859-1'))


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Merge analysis results from parallel analysis.')

    parser.add_argument('--working_dir', dest='wdir', type=str, required=True,
                        help = "The name of the analysis working directory.")
    parser.add_argument('--bin', dest='merged', type=str, required=True,
                        help = "The name of the merged localization file.")
    parser.add_argument('--ext', dest='ext', type=str, required=False, default=".bin",
                        help = "The name of an extension, if any.")

    args = parser.parse_args()

    mergeAnalysis(args.wdir, args.merged, args.ext)
