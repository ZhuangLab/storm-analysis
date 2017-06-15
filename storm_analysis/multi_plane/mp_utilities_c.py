#!/usr/bin/env python
"""
Utility functions for multi-plane fitting. 

Hazen 06/17
"""
import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import pickle

import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.loadclib as loadclib


# C library loading & function definitions.

mp_util = loadclib.loadCLibrary("storm_analysis.multi_plane", "mp_utilities")

mp_util.mpuBadPeakMask.argtypes = [ctypes.c_void_p,
                                   ndpointer(dtype=numpy.float64),
                                   ndpointer(dtype=numpy.uint8),
                                   ctypes.c_int]

mp_util.mpuCleanup.argtypes = [ctypes.c_void_p]

mp_util.mpuFilterPeaks.argtypes = [ctypes.c_void_p,
                                   ndpointer(dtype=numpy.float64),
                                   ndpointer(dtype=numpy.float64),
                                   ndpointer(dtype=numpy.uint8),
                                   ctypes.c_int,
                                   ctypes.c_int]

mp_util.mpuInitialize.argtypes = [ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_int,
                                  ctypes.c_int,
                                  ctypes.c_int,
                                  ctypes.c_int]
                                  
mp_util.mpuInitialize.restype = ctypes.c_void_p

mp_util.mpuMarkClosePeaks.argtypes = [ctypes.c_void_p,
                                      ndpointer(dtype=numpy.float64),
                                      ndpointer(dtype=numpy.uint8),
                                      ctypes.c_int]

mp_util.mpuMergeNewPeaks.argtypes = [ctypes.c_void_p,
                                     ndpointer(dtype=numpy.float64),
                                     ndpointer(dtype=numpy.float64),
                                     ndpointer(dtype=numpy.float64),
                                     ctypes.c_int,
                                     ctypes.c_int]

mp_util.mpuSetTransforms.argtypes = [ctypes.c_void_p,
                                     ndpointer(dtype=numpy.float64),
                                     ndpointer(dtype=numpy.float64)]
                                     
mp_util.mpuSplitPeaks.argtypes = [ctypes.c_void_p,
                                  ndpointer(dtype=numpy.float64),
                                  ndpointer(dtype=numpy.float64),
                                  ctypes.c_int]


class MpUtil(object):

    def __init__(self,
                 radius = None,
                 neighborhood = None,
                 im_size_x = None,
                 im_size_y = None,
                 n_channels = None,
                 n_zplanes = None,
                 **kwds):
        super().__init__(**kwds)
        self.n_channels = n_channels
        self.mpu = mp_util.mpuInitialize(radius,
                                         neighborhood,
                                         im_size_x,
                                         im_size_y,
                                         n_channels,
                                         n_zplanes)

    def badPeakMask(self, peaks):
        """
        Return a mask for that can be passed to filterPeaks() to 
        remove bad (ERROR) peaks.
        """
        mask = numpy.ascontiguousarray(numpy.ones(peaks.shape[0], dtype = numpy.uint8))
        mp_util.mpuBadPeakMask(self.mpu,
                               numpy.ascontiguousarray(peaks, dtype = numpy.float64),
                               mask,
                               mask.size)
        return mask

    def cleanup(self):
        mp_util.mpuCleanup(self.mpu)

    def filterPeaks(self, peaks, mask):
        assert((peaks.shape[0] % self.n_channels) == 0)
        assert((mask.size * self.n_channels) == peaks.shape[0])

        in_peaks_size = int(peaks.shape[0]/self.n_channels)
        out_peaks_size = numpy.count_nonzero(mask)
        out_peaks = numpy.ascontiguousarray(numpy.zeros((out_peaks_size*self.n_channels, utilC.getNPeakPar())),
                                            dtype = numpy.float64)

        mp_util.mpuFilterPeaks(self.mpu,
                               numpy.ascontiguousarray(peaks, dtype = numpy.float64),
                               out_peaks,
                               numpy.ascontiguousarray(mask, dtype = numpy.uint8),
                               in_peaks_size,
                               out_peaks_size)
        return out_peaks

    def markClosePeaks(self, peaks):
        assert((peaks.shape[0] % self.n_channels) == 0)

        mask_size = int(peaks.shape[0]/self.n_channels)
        peaks = numpy.ascontiguousarray(peaks, dtype = numpy.float64)
        mask = numpy.ascontiguousarray(numpy.ones(mask_size, dtype = numpy.uint8))
        
        mp_util.mpuMarkClosePeaks(self.mpu,
                                  peaks,
                                  mask,
                                  mask_size)
        return [peaks, mask]
        
    def mergeNewPeaks(self, cur_peaks, new_peaks):
        assert((cur_peaks.shape[0] % self.n_channels) == 0)
        assert((new_peaks.shape[0] % self.n_channels) == 0)

        n_cur = int(cur_peaks.shape[0]/self.n_channels)
        n_new = int(new_peaks.shape[0]/self.n_channels)
        merged_peaks = numpy.ascontiguousarray(numpy.zeros(((n_cur + n_new)*self.n_channels, utilC.getNPeakPar())),
                                               dtype = numpy.float64)

        # Create merged peaks list.
        mp_util.mpuMergeNewPeaks(self.mpu,
                                 numpy.ascontiguousarray(cur_peaks, dtype = numpy.float64),
                                 numpy.ascontiguousarray(new_peaks, dtype = numpy.float64),
                                 merged_peaks,
                                 n_cur,
                                 n_new)

        # Create mask for bad (i.e. ERROR) peaks.
        mask = self.badPeakMask(merged_peaks)

        # Return filtered merged_peaks.
        return self.filterPeaks(merged_peaks, mask)

    def prettyPrintPeak(self, peaks, index):
        assert((peaks.shape[0] % self.n_channels) == 0)
    
        i_b = utilC.getBackgroundIndex()
        i_h = utilC.getHeightIndex()
        i_x = utilC.getXCenterIndex()
        i_y = utilC.getYCenterIndex()
        i_z = utilC.getZCenterIndex()

        n_peaks = int(peaks.shape[0]/self.n_channels)
        for i in range(self.n_channels):
            j = i*n_peaks + index
            print("{0:d}) {1:.2f} {2:.2f} {3:.2f} {4:.2f} {5:.2f}".format(i,
                                                                          peaks[j,i_b],
                                                                          peaks[j,i_h],
                                                                          peaks[j,i_x],
                                                                          peaks[j,i_y],
                                                                          peaks[j,i_z]))

    def removeClosePeaks(self, peaks):
        """
        If two peaks are closer than radius, remove the dimmer 
        one and mark the brighter one as running.
        """
        return self.filterPeaks(*self.markClosePeaks(peaks))

    def setTransforms(self, xt, yt):
        """
        Set channel0 to channelN transforms.
        """
        mp_util.mpuSetTransforms(self.mpu,
                                 numpy.ascontiguousarray(xt, dtype = numpy.float64),
                                 numpy.ascontiguousarray(yt, dtype = numpy.float64))

    def splitPeaks(self, peaks):
        """
        Create per-channel peaks from channel0 peaks. Note that you
        must set the channel to channel transforms first with setTransforms().
        """
        split_peaks = numpy.ascontiguousarray(numpy.zeros((peaks.shape[0]*self.n_channels, utilC.getNPeakPar())),
                                              dtype = numpy.float64)
        mp_util.mpuSplitPeaks(self.mpu,
                              numpy.ascontiguousarray(peaks, dtype = numpy.float64),
                              split_peaks,
                              peaks.shape[0])
        return split_peaks

    def testCreatePeaks(self, x, y, z = None):
        """
        Create peak arrays for testing purposes.
        """
        assert(x.size == y.size)
        if z is not None:
            assert(x.size == z.size)

        i_s = utilC.getStatusIndex()
        i_x = utilC.getXCenterIndex()
        i_y = utilC.getYCenterIndex()
        i_z = utilC.getZCenterIndex()

        peaks = numpy.zeros((x.size, utilC.getNPeakPar()))
        for i in range(x.size):
            peaks[i,:] = float(i+1)
            peaks[i,i_s] = 0.0
            peaks[i,i_x] = x[i]
            peaks[i,i_y] = y[i]
            if z is not None:
                peaks[i,i_z] = z[i]

        return peaks
                

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
        #peaks[i,i_z] = 0.01 * (i%4)
        #peaks[i,i_z] = zi
        peaks[i,i_z] = z_values[zi]

def loadMappings(filename, margin):

    max_ch = 1
    mappings = {}

    if filename is not None:
        # Load the mappings hash.
        with open(filename, 'rb') as fp:
            mappings = pickle.load(fp)

        # Figure out number of channels.
        while "0_" + str(max_ch) + "_x" in mappings:
            max_ch += 1

    # Create and fill in mapping transform arrays.
    xt_0toN = numpy.zeros((max_ch,3))
    yt_0toN = numpy.zeros((max_ch,3))
    xt_Nto0 = numpy.zeros((max_ch,3))
    yt_Nto0 = numpy.zeros((max_ch,3))

    # 0 <-> 0 is the identity transform.
    xt_0toN[0,1] = 1.0
    yt_0toN[0,2] = 1.0
    xt_Nto0[0,1] = 1.0
    yt_Nto0[0,2] = 1.0

    for i in range(1,max_ch):
        xt_0toN[i,:] = marginCorrect(mappings["0_" + str(i) + "_x"], margin)
        yt_0toN[i,:] = marginCorrect(mappings["0_" + str(i) + "_y"], margin)
        xt_Nto0[i,:] = marginCorrect(mappings[str(i) + "_0_x"], margin)
        yt_Nto0[i,:] = marginCorrect(mappings[str(i) + "_0_y"], margin)

    xt_0toN = numpy.ascontiguousarray(xt_0toN, dtype = numpy.float64)
    yt_0toN = numpy.ascontiguousarray(yt_0toN, dtype = numpy.float64)
    xt_Nto0 = numpy.ascontiguousarray(xt_Nto0, dtype = numpy.float64)
    yt_Nto0 = numpy.ascontiguousarray(yt_Nto0, dtype = numpy.float64)
    
    return [xt_0toN, yt_0toN, xt_Nto0, yt_Nto0]
    
def marginCorrect(tr, margin):
    """
    Correct affine transform for the margin that was added to an image.
    """
    tr[0] += margin - (tr[1] + tr[2]) * margin
    return tr



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
