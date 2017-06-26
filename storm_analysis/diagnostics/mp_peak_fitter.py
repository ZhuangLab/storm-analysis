#!/usr/bin/env python
"""
Evaluate the performance of the multi-plane peak fitter. How
good is it? What is the optimal way to weight the different
planes as a function of Z?

Designed to be used in combination with sa_utilities.finding_fitting_error.py

Hazen 06/17
"""
import numpy

import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.parameters as params
import storm_analysis.sa_library.readinsight3 as readinsight3

import storm_analysis.sa_utilities.std_analysis as stdAnalysis

import storm_analysis.multi_plane.find_peaks_std as findPeaksStd


xy_error = None
z_error = None

class PFMovieReader(findPeaksStd.MPMovieReader):

    def __init__(self, bin_name = None, **kwds):
        super().__init__(**kwds)

        self.i3_data = readinsight3.loadI3File(bin_name)

    def getLocalizations(self):
        
        #
        # Return groud truth localization positions in the multi-fit
        # peak format.
        #
        mf_peaks = i3dtype.convertToMultiFit(self.i3_data,
                                             self.movie_x,
                                             self.movie_y,
                                             self.cur_frame,
                                             self.nm_per_pixel)
        return mf_peaks
        
    
class PFPeakFinder(findPeaksStd.MPPeakFinder):

    def findPeaks(self, peaks):

        # Adjust x,y for margin.
        peaks[:,utilC.getXCenterIndex()] += float(self.margin)
        peaks[:,utilC.getYCenterIndex()] += float(self.margin)
        
        # Convert z (in nano-meters) to spline z units.
        zi = utilC.getZCenterIndex()
        peaks[:,zi] = self.s_to_psfs[0].getScaledZ(peaks[:,zi])

        # Split peaks into per-channel peaks.
        return self.mpu.splitPeaks(peaks)
    
    
class PFPeakFinderFitter(findPeakStd.MPPeakFinderFitter):

   def analyzeImage(self, movie_reader, save_residual = False, verbose = False):

       # Load images & background estimates.
       [images, fit_peaks_images] = self.loadImages(movie_reader)

       # Pass images to fitter.
       self.peak_fitter.newImages(images)

       # Load ground truth localization positions.
       truth_peaks = movieReader.getLocalizations()

       # Adjust to work with multi-plane.
       fit_peaks = self.peak_finder.findPeaks(truth_peaks.copy())

       # Calculate best fits peaks.
       [fit_peaks, fit_peaks_images] = self.peak_fitter.fitPeaks(fit_peaks)

       # Adjust for margin.
       fit_peaks[:,utilC.getXCenterIndex()] -= float(self.margin)
       fit_peaks[:,utilC.getYCenterIndex()] -= float(self.margin)
       
       return [fit_peaks, None]
        
