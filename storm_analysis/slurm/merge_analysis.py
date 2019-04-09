#!/usr/bin/env python
"""
Merges the intermediate localization files into a single 
localization file.

Hazen 08/17
"""

import glob
import os
from xml.etree import ElementTree

import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.slurm.check_analysis as checkAnalysis

def mergeAnalysis(dir_name, h5_name):

    with saH5Py.SAH5Py(h5_name, is_existing = False, sa_type = "merged") as merged_h5:

        # Get job XML files.
        job_xml_files = checkAnalysis.getSortedJobXML(dir_name)
        
        # Check for corresponding HDF5 files.
        job_complete = True
        for i in range(len(job_xml_files)):

            sub_h5_name = os.path.join(dir_name, "p_" + str(i+1) + ".hdf5")
            
            if os.path.exists(sub_h5_name):
                with saH5Py.SAH5Py(sub_h5_name) as h5:
                    if not h5.isAnalysisFinished():
                        job_complete = False
                        break
                    
                    if (i == 0):
                        merged_h5.setMovieInformation(*h5.getMovieInformation())
                        merged_h5.setPixelSize(h5.getPixelSize())
                        
                        # Use XML metadata from the first file.
                        merged_h5.addMetadata(h5.getMetadata())

                    # Copy localizations. There shouldn't be any tracking
                    # information at this stage of the analysis.
                    for fnum, locs in h5.localizationsIterator():
                        merged_h5.addLocalizations(locs, fnum)

            else:
                job_complete = False
                break

        if not job_complete:
            print("Merge failed because", job_xml_files[i], "is incomplete.")


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Merge analysis results from parallel analysis.')

    parser.add_argument('--working_dir', dest='wdir', type=str, required=True,
                        help = "The name of the analysis working directory.")
    parser.add_argument('--h5_name', dest='merged', type=str, required=True,
                        help = "The name for the merged localization file.")
    parser.add_argument('--ext', dest='ext', type=str, required=False, default=[".bin"], nargs = "*",
                        help = "The name of the extensions, if any.")

    args = parser.parse_args()

    mergeAnalysis(args.wdir, args.merged, args.ext)
