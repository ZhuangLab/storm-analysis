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

import storm_analysis.sa_library.sa_h5py as saH5Py


def checkAnalysis(dir_name):

    job_xml_files = getSortedJobXML(dir_name)

    # Check for corresponding HDF5 files.
    incomplete = None
    for i in range(len(job_xml_files)):

        if ((i%20)==0):
            print("Checking", job_xml_files[i])

        h5_name = os.path.join(dir_name, "p_" + str(i+1) + ".hdf5")
        if os.path.exists(h5_name):
            with saH5Py.SAH5Py(h5_name) as h5:
                if h5.isAnalysisFinished():
                    continue

        print("Job", job_xml_files[i], "is incomplete.")
        if incomplete is None:
            incomplete = str(i+1)
        else:
            incomplete += "," + str(i+1)

    if incomplete is not None:
        print("suggested job array string:")
        print(incomplete)


def getSortedJobXML(dir_name):
    
    # Find all the job*.xml files.
    job_xml_files = glob.glob(os.path.join(dir_name,"job*.xml"))

    # Sort job files.
    job_xml_files = sorted(job_xml_files, key = lambda x: int(os.path.splitext(os.path.basename(x))[0].split("_")[1]))

    return job_xml_files


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Check analysis results for parallel analysis.')

    parser.add_argument('--working_dir', dest='wdir', type=str, required=True,
                        help = "The name of the analysis working directory.")

    args = parser.parse_args()

    checkAnalysis(args.wdir)

