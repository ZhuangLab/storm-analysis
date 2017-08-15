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

def mergeAnalysis(dir_name, bin_base_name, extensions = [".bin"]):

    # Create Insight3 file writers.
    i3_out = []
    for ext in extensions:
        i3_out.append(writeinsight3.I3Writer(bin_base_name + ext))

    # Find all the job*.xml files.
    job_xml_files = glob.glob(dir_name + "job*.xml")

    # Sort job files.
    job_xml_files = sorted(job_xml_files, key = lambda x: int(os.path.splitext(os.path.basename(x))[0].split("_")[1]))

    # Check for corresponding mlist.bin files.
    metadata = None
    last_frame = 0
    for i in range(len(job_xml_files)):

        job_complete = True
        for j, ext in enumerate(extensions):
            mlist_name = dir_name + "p_" + str(i+1) + "_mlist" + ext

            if os.path.exists(mlist_name) and readinsight3.checkStatus(mlist_name):

                # Load metadata from the first file.
                if (i == 0) and (j == 0):
                    metadata = readinsight3.loadI3Metadata(mlist_name)

                # Read localizations.
                i3_data = readinsight3.loadI3File(mlist_name, verbose = False)

                # Print frame range covered.
                if (j == 0):
                    last_frame = i3_data["fr"][-1]
                    print(i3_data["fr"][0], last_frame, mlist_name)

                # Add localizations to the output file.
                i3_out[j].addMolecules(i3_data)

            else:
                job_complete = False
                break


        if not job_complete:
            print("Merge failed because", job_xml_files[i], "is incomplete.")
            for j, ext in enumerate(extensions):
                i3_out[j].close()
                os.remove(bin_base_name + ext)
            assert(False)

    if metadata is None:
        print("No metadata found.")
        for i3w in i3_out:
            i3w.close()
    else:

        # Fix movie length node based on the last frame of the last molecule.
        metadata.find("movie").find("movie_l").text = str(last_frame)

        # Also need to fix analysis end points. We are assuming that the
        # entire movie was analyzed.
        metadata.find("settings").find("start_frame").text = "-1"
        metadata.find("settings").find("max_frame").text = "-1"

        for i3w in i3_out:
            i3w.closeWithMetadata(ElementTree.tostring(metadata, 'ISO-8859-1'))


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Merge analysis results from parallel analysis.')

    parser.add_argument('--working_dir', dest='wdir', type=str, required=True,
                        help = "The name of the analysis working directory.")
    parser.add_argument('--bin_base_name', dest='merged', type=str, required=True,
                        help = "The base name of the merged localization file (i.e. without .bin extension)")
    parser.add_argument('--ext', dest='ext', type=str, required=False, default=[".bin"], nargs = "*",
                        help = "The name of the extensions, if any.")

    args = parser.parse_args()

    mergeAnalysis(args.wdir, args.merged, args.ext)
