#!/usr/bin/env python
"""
Creates an analysis.xml for file for a specific range of frames. This
is for the purpose of being able to run analysis in parallel on a 
single STORM movie.

Hazen 07/17
"""
import os

from xml.etree import ElementTree


def splitAnalysisXML(working_dir, params_xml, start_frame, max_frame, divisions):
    """
    working_dir - The working directory where the analysis will be performed.
    params_xml - The name of a parameters XML file.
    start_frame - Which frame to start the analysis on.
    max_frame - Which frame to stop the analysis on.
    divisions - How many chunks to break the movie into.
    """
    assert(divisions >= 2)
    
    # Load the original analysis XML file.
    params = ElementTree.parse(params_xml).getroot()

    # Get relevant nodes.
    start_frame_node = params.find("start_frame")
    max_frame_node = params.find("max_frame")
        
    # Set radius to 0.0 to disable tracking.
    params.find("radius").text = "0.0"

    # Turn off drift correction.
    params.find("drift_correction").text = "0"

    # Create new analysis files for each division.
    step_size = int(max_frame/divisions) + 1

    index = 1
    cur_frame = start_frame
    while (cur_frame < max_frame):
        
        # The first frames can have lots of localizations so we don't
        # want this job to have a lot of frames in it.
        if (index == 1) and (step_size > 10):
            stop_frame = 10
        else:
            stop_frame = cur_frame + step_size

        # Don't go past the maximum frame.
        if (stop_frame >= max_frame):
            stop_frame = max_frame

        # Some feeback.
        if ((index % 20) == 0):
            print("Creating XML for job", index, start_frame, stop_frame)

        start_frame_node.text = str(cur_frame)
        max_frame_node.text = str(stop_frame)
            
        # Save XML.
        with open(os.path.join(working_dir, "job_" + str(index) + ".xml"), "wb") as fp:
            fp.write(ElementTree.tostring(params, 'ISO-8859-1'))

        index += 1
        if (stop_frame == max_frame):
            cur_frame = max_frame + 1
        else:
            cur_frame = stop_frame

    print("Created", index - 1, "jobs.")
        
    
if (__name__ == "__main__"):

    import argparse

    import storm_analysis.sa_library.datareader as datareader

    parser = argparse.ArgumentParser(description = 'Split analysis XML for parallel analysis.')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie that will be analyzed.")
    parser.add_argument('--working_dir', dest='wdir', type=str, required=True,
                        help = "The directory for intermediate analysis.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")
    parser.add_argument('--divisions', dest='divisions', type=int, required=True,
                        help = "How many sections to break the movie into.")

    args = parser.parse_args()

    # Figure out how many frames are in the movie.
    movie = datareader.inferReader(args.movie)
    movie_length = movie.filmSize()[2]

    print("Movie has", movie_length, "frames")
    splitAnalysisXML(args.wdir, args.settings, 0, movie_length, args.divisions)


