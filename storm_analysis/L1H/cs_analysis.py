#!/usr/bin/env python
"""
Perform compressed sensing analysis on a dax file using the
homotopy approach. Return the results in hres image format and
as a list of object locations.

Hazen 09/12
"""

import numpy

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.parameters as parameters
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

import storm_analysis.L1H.setup_A_matrix as setup_A_matrix
import storm_analysis.L1H.homotopy_imagea_c as homotopy_imagea_c


def analyze(movie_name, settings_name, hres_name, bin_name):

    movie_data = datareader.inferReader(movie_name)

    #
    # FIXME:
    #
    # This should also start at the same frame as hres in the event of a restart.
    #
    i3_file = writeinsight3.I3Writer(bin_name)
    
    params = parameters.ParametersL1H().initFromFile(settings_name)

    #
    # Load the a matrix and setup the homotopy image analysis class.
    #
    a_mat_file = params.getAttr("a_matrix")

    print("Using A matrix file:", a_mat_file)
    a_mat = setup_A_matrix.loadAMatrix(a_mat_file)

    image = movie_data.loadAFrame(0)
    htia = homotopy_imagea_c.HomotopyIA(a_mat,
                                        params.getAttr("epsilon"),
                                        image.shape)

    #
    # This opens the file. If it already exists, then it sets the file pointer
    # to the end of the file & returns the number of the last frame analyzed.
    #
    curf = htia.openHRDataFile(hres_name)

    #
    # Figure out which frame to start & stop at.
    #
    [dax_x,dax_y,dax_l] = movie_data.filmSize()

    if params.hasAttr("start_frame"):
        if (params.getAttr("start_frame") >= curf) and (params.getAttr("start_frame") < dax_l):
            curf = params.getAttr("start_frame")

    if params.hasAttr("max_frame"):
        if (params.getAttr("max_frame") > 0) and (params.getAttr("max_frame") < dax_l):
            dax_l = params.getAttr("max_frame")

    print("Starting analysis at frame", curf)

    #
    # Analyze the dax data.
    #
    total_peaks = 0
    try:
        while(curf<dax_l):

            # Load image, subtract baseline & remove negative values.
            image = movie_data.loadAFrame(curf).astype(numpy.float)

            # Convert to photo-electrons.
            image -= params.getAttr("camera_offset")
            image = image/params.getAttr("camera_gain")

            # Remove negative values.
            mask = (image < 0)
            image[mask] = 0

            # Analyze image.
            hres_image = htia.analyzeImage(image)
            peaks = htia.saveHRFrame(hres_image, curf + 1)
            [cs_x,cs_y,cs_a,cs_i] = htia.getPeaks(hres_image)
            i3_file.addMoleculesWithXYAItersFrame(cs_x, cs_y, cs_a, cs_i, curf+1)

            peaks = cs_x.size
            total_peaks += peaks
            print("Frame:", curf, peaks, total_peaks)

            curf += 1

    except KeyboardInterrupt:
        print("Analysis stopped.")

    # cleanup
    htia.closeHRDataFile()
    i3_file.close()


if (__name__ == "__main__"):
    
    import argparse

    parser = argparse.ArgumentParser(description = 'L1H analysis - Babcock, Optics Express, 2013')

    parser.add_argument('--movie', dest='movie', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")
    parser.add_argument('--hres', dest='hres', type=str, required=True,
                        help = "The name of 'high resolution' output file. This a compressed version of the final image.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file. This is a binary file in Insight3 format.")

    args = parser.parse_args()
    
    analyze(args.movie, args.settings, args.hres, args.mlist)
    
#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
