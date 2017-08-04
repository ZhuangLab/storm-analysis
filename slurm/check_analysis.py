#!/usr/bin/env python
"""
Checks that the analysis is complete. Basically this looks
at the all analysis.xml files and checks that there is a
corresponding mlist.bin file and that this file is not corrupt
as verified based on the status byte.

Hazen 08/17
"""

import glob
import os

import storm_analysis.sa_library.readinsight3 as readinsight3

def checkAnalysis(dir_name):
    
    # Find all the job*.xml files.
    job_xml_files = glob.glob(dir_name + "job*.xml")

    # Sort job files.
    job_xml_files = sorted(job_xml_files, key = lambda x: int(os.path.splitext(os.path.basename(x))[0].split("_")[1]))

    # Check for corresponding mlist.bin files.
    for i in range(len(job_xml_files)):

        mlist_name = dir_name + "p_" + str(i+1) + "_mlist.bin"
        if os.path.exists(mlist_name):
            if readinsight3.checkStatus(mlist_name):
                continue

        print("Job", job_xml_files[i], "is incomplete.")


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Check analysis results for parallel analysis.')

    parser.add_argument('--working_dir', dest='wdir', type=str, required=True,
                        help = "The name of the analysis working directory.")

    args = parser.parse_args()

    checkAnalysis(args.wdir)

