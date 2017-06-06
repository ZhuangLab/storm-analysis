#!/usr/bin/env python
"""
Utility functions for multi-plane fitting. 

This does not actually interface with a C library yet, but 
likely will soon once we start trying to optimize this.

Hazen 06/17
"""
import numpy

import storm_analysis.sa_library.ia_utilities_c as utilC


def initializeBackground(peaks, backgrounds):
    """
    Set initial background values for the peaks.
    """
    i_b = utilC.getBackgroundIndex()
    i_x = utilC.getXCenterIndex()
    i_y = utilC.getYCenterIndex()

    n_channels = len(backgrounds)
    assert((peaks.shape[0] % n_channels) == 0)

    n_peaks = int(peaks.shape[0]/n_channels)
    for i in range(n_peaks):
        for j in range(n_channels):
            k = i*n_channels + j
            xi = int(round(peaks[k,i_x]))
            yi = int(round(peaks[k,i_y]))
            peaks[k,i_b] = backgrounds[j][yi,xi]

def initializeHeight(peaks, foregrounds, height_rescale):
    """
    Set initial height values for the peaks.
    """
    i_h = utilC.getHeightIndex()
    i_x = utilC.getXCenterIndex()
    i_y = utilC.getYCenterIndex()
    i_z = utilC.getZCenterIndex()

    n_channels = len(foregrounds[0])
    assert((peaks.shape[0] % n_channels) == 0)

    n_peaks = int(peaks.shape[0]/n_channels)
    for i in range(n_peaks):
        for j in range(n_channels):
            k = i*n_channels + j
            xi = int(round(peaks[k,i_x]))
            yi = int(round(peaks[k,i_y]))
            zi = int(round(peaks[k,i_z]))
            
            peaks[k,i_h] = foregrounds[zi][j][yi,xi] * height_rescale[zi][j]

def initializeZ(peaks, z_values):
    """
    Set initial Z values for the peaks.
    """
    i_z = utilC.getZCenterIndex()
    
    for i in range(peaks.shape[0]):
        zi = int(round(peaks[i,i_z]))
        peaks[i,i_z] = 0.01 * (i%4)
        #peaks[i,i_z] = zi
        #peaks[i,i_z] = z_values[zi]

def marginCorrect(tr, margin):
    """
    Correct affine transform for the margin that was added to an image.
    """
    tr[0] += margin - (tr[1] + tr[2]) * margin
    return tr

def mergeNewPeaks(peak, new_peaks, radius, neighborhood, n_channels):
    assert(False)
    return new_peaks

def prettyPrintPeak(peaks, index, n_channels):
    
    i_b = utilC.getBackgroundIndex()
    i_h = utilC.getHeightIndex()
    i_x = utilC.getXCenterIndex()
    i_y = utilC.getYCenterIndex()
    i_z = utilC.getZCenterIndex()
    
    for i in range(n_channels):
        j = index*n_channels + i
        print("{0:d}) {1:.2f} {2:.2f} {3:.2f} {4:.2f} {5:.2f}".format(i,
                                                                      peaks[j,i_b],
                                                                      peaks[j,i_h],
                                                                      peaks[j,i_x],
                                                                      peaks[j,i_y],
                                                                      peaks[j,i_z]))

def splitPeaks(peaks, xt, yt):
    """
    Split peaks array into channels and map x and y locations.
    """
    n_channels = len(xt) + 1
    i_x = utilC.getXCenterIndex()
    i_y = utilC.getYCenterIndex()

    out_peaks = numpy.zeros((peaks.shape[0]*n_channels, peaks.shape[1]))
    for i in range(peaks.shape[0]):
        for j in range(n_channels):
            k = i*n_channels

            # First just copy all the values.
            out_peaks[k+j,:] = peaks[i,:].copy()

            if (j > 0):
                xtj = xt[j-1]
                ytj = yt[j-1]

                #
                # Then correct x, y positions.
                #
                # Need to transpose x and y because the mapping was
                # determined using the transposes of the channel images.
                #
                yi = peaks[i,i_x]
                xi = peaks[i,i_y]
                yf = xtj[0] + xtj[1]*xi + xtj[2]*yi
                xf = ytj[0] + ytj[1]*xi + ytj[2]*yi

                out_peaks[k+j,i_x] = xf
                out_peaks[k+j,i_y] = yf

    # This should break the connection between the input and
    # the output peaks by forcing a copy operation.
    return numpy.ascontiguousarray(out_peaks)
