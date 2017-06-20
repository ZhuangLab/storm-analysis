#!/usr/bin/env python
"""
Evaluate the performance of the multi-plane peak finder.

Specifically, how close are the guesses that it provides to 
the actual values after fitting?

Hazen 06/17
"""
import numpy

import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_utilities.std_analysis as stdAnalysis

import storm_analysis.multi_plane.find_peaks_std as findPeaksStd

fitted_peaks = None
found_peaks = None

class PFPeakFinder(findPeaksStd.MPPeakFinder):

    def findPeaks(self, no_bg_image, peaks):
        [found_new_peaks, peaks] = super().findPeaks(no_bg_image, peaks)
        
        global found_peaks
        found_peaks = peaks.copy()
        
        return [found_new_peaks, peaks]
        

class PFPeakFitter(findPeaksStd.MPPeakFitter):

    def fitPeaks(self, peaks):
        [fit_peaks, fit_peaks_images] = super().fitPeaks(peaks)

        global fitted_peaks
        fitted_peaks = fit_peaks.copy()
        
        return [fit_peaks, fit_peaks_images]

        
def analyze(base_name, mlist_name, settings_name):
    parameters = params.ParametersMultiplane().initFromFile(settings_name)

    # Set to only analyze 1 frame, 1 iteration of peak finding.
    parameters.setAttr("max_frame", "int", 1)
    parameters.setAttr("iterations", "int", 1)

    # Analyze one frame.
    finder = PFPeakFinder(parameters)
    fitter = PFPeakFitter(parameters)
    finder_fitter = findPeaksStd.MPFinderFitter(parameters, finder, fitter)
    reader = findPeaksStd.MPMovieReader(base_name = base_name,
                                        parameters = parameters)
    stdAnalysis.standardAnalysis(finder_fitter,
                                 reader,
                                 mlist_name,
                                 parameters)

    # Evaluate found vs fit parameters.
    hi = utilC.getHeightIndex()
    xi = utilC.getXCenterIndex()
    yi = utilC.getYCenterIndex()

    # First identify peak pairs.
    fi_x = fitted_peaks[:,xi]
    fi_y = fitted_peaks[:,yi]
    fo_x = found_peaks[:,xi]
    fo_y = found_peaks[:,yi]

    dd = utilC.peakToPeakDist(fo_x, fo_y, fi_x, fi_y)
    di = utilC.peakToPeakIndex(fo_x, fo_y, fi_x, fi_y)
    
    print(numpy.mean(dd), fi_x.size, fo_x.size)

    fo_h = []
    fi_h = []
    for i in range(dd.size):
        if(dd[i]<2.0):
            fo_h.append(found_peaks[i,hi])
            fi_h.append(fitted_peaks[di[i],hi])

    fi_h = numpy.array(fi_h)
    fo_h = numpy.array(fo_h)

    print(numpy.mean(fi_h), fi_h.size)
    print(numpy.mean(fo_h), fo_h.size)
    print(numpy.mean(fi_h)/numpy.mean(fo_h))
    

if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Peak finder diagnostics - ...')

    parser.add_argument('--basename', dest='basename', type=str, required=True,
                        help = "The base name of the movie to analyze. Movies can be in .dax, .tiff or .spe format.")
    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations output file. This is a binary file in Insight3 format.")
    parser.add_argument('--xml', dest='settings', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()
    
    analyze(args.basename, args.mlist, args.settings)
