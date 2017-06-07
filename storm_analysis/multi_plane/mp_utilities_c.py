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
    for i in range(n_channels):
        for j in range(n_peaks):
            k = i*n_peaks + j
            xi = int(round(peaks[k,i_x]))
            yi = int(round(peaks[k,i_y]))
            peaks[k,i_b] = backgrounds[i][yi,xi]

def initializeHeight(peaks, foregrounds, height_rescale):
    """
    Set initial height values for the peaks.
    """
    i_h = utilC.getHeightIndex()
    i_x = utilC.getXCenterIndex()
    i_y = utilC.getYCenterIndex()
    i_z = utilC.getZCenterIndex()

    n_channels = len(foregrounds[0])
    inv_n_channels = 1.0/n_channels
    assert((peaks.shape[0] % n_channels) == 0)

    n_peaks = int(peaks.shape[0]/n_channels)
    for i in range(n_peaks):

        #
        # Compute average height.
        #
        height = 0.0
        for j in range(n_channels):
            k = j*n_peaks + i
            xi = int(round(peaks[k,i_x]))
            yi = int(round(peaks[k,i_y]))
            zi = int(round(peaks[k,i_z]))
            height += foregrounds[zi][j][yi,xi] * height_rescale[zi][j]

        #
        # Assign the same height to every peak in the group. Throughout
        # the analysis we try and make sure that all the peaks in a group
        # have the same height. We are assuming that the splines were
        # properly normalized so that this makes sense.
        #
        for j in range(n_channels):
            k = j*n_peaks + i
            peaks[k,i_h] = height * inv_n_channels

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

    n_peaks = int(peaks.shape[0]/n_channels)
    for i in range(n_channels):
        j = i*n_peaks + index
        print("{0:d}) {1:.2f} {2:.2f} {3:.2f} {4:.2f} {5:.2f}".format(i,
                                                                      peaks[j,i_b],
                                                                      peaks[j,i_h],
                                                                      peaks[j,i_x],
                                                                      peaks[j,i_y],
                                                                      peaks[j,i_z]))

def removeClosePeaks(peaks, sigma, neighborhood):
    return peaks

def splitPeaks(peaks, xt, yt):
    """
    Split peaks array into channels and map x and y locations.
    """
    n_channels = len(xt) + 1
    i_x = utilC.getXCenterIndex()
    i_y = utilC.getYCenterIndex()

    n_peaks = peaks.shape[0]
    out_peaks = numpy.zeros((n_peaks*n_channels, peaks.shape[1]))
    for i in range(n_channels):
        for j in range(n_peaks):
            k = i*n_peaks

            # First just copy all the values.
            out_peaks[k+j,:] = peaks[j,:].copy()

            #
            # Then correct x, y positions for peaks that are not in
            # channel 0.
            #
            # Need to transpose x and y because the mapping was
            # determined using the transposes of the channel images.
            #            
            if (i > 0):
                xtj = xt[i-1]
                ytj = yt[i-1]

                yi = peaks[j,i_x]
                xi = peaks[j,i_y]
                yf = xtj[0] + xtj[1]*xi + xtj[2]*yi
                xf = ytj[0] + ytj[1]*xi + ytj[2]*yi

                out_peaks[k+j,i_x] = xf
                out_peaks[k+j,i_y] = yf

    # This should break the connection between the input and
    # the output peaks by forcing a copy operation.
    return numpy.ascontiguousarray(out_peaks)
