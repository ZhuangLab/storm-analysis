#!/usr/bin/env python
"""
Tests for multi_plane.mp_utilities_c
"""

import numpy
import pickle

import storm_analysis

import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.multi_plane.mp_utilities_c as mpUtilC

def createMPU(n_channels):
    return mpUtilC.MpUtil(radius = 1.0,
                          neighborhood = 5.0,
                          im_size_x = 200,
                          im_size_y = 100,
                          n_channels = n_channels,
                          n_zplanes = 1)

def createTransforms(n_channels):
    xt = numpy.zeros((n_channels, 3))
    yt = numpy.zeros((n_channels, 3))
    xt[:,1] = 1.0
    yt[:,2] = 1.0
    if(n_channels>1):
        xt[1:,1] = -1.0
        yt[1:,2] = -1.0
    return [xt,yt]


def test_bad_peak_mask_1():
    mpu = createMPU(1)
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                numpy.array([10.0, 10.0]))
    mask = mpu.badPeakMask(peaks)
    assert(mask[0] == 1)

def test_bad_peak_mask_2():
    mpu = createMPU(1)
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                numpy.array([10.0, 10.0]))
    i_s = utilC.getStatusIndex()
    peaks[:,i_s] = 2.0

    mask = mpu.badPeakMask(peaks)
    assert(mask[0] == 0)

    
def test_filter_1():
    mpu = createMPU(1)
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                numpy.array([10.0, 10.0]))
    mask = numpy.array([True, True])
    fpeaks = mpu.filterPeaks(peaks, mask)

    assert(fpeaks.shape[0] == 2)

def test_filter_2():
    mpu = createMPU(1)
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                numpy.array([10.0, 10.0]))
    mask = numpy.array([True, False])
    fpeaks = mpu.filterPeaks(peaks, mask)

    assert(fpeaks.shape[0] == 1)
    assert(fpeaks[0,0] == 1.0)

def test_filter_3():
    mpu = createMPU(1)
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                numpy.array([10.0, 10.0]))
    mask = numpy.array([False, True])
    fpeaks = mpu.filterPeaks(peaks, mask)
    
    assert(fpeaks.shape[0] == 1)
    assert(fpeaks[0,0] == 2.0)

def test_filter_4():
    mpu = createMPU(1)
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                numpy.array([10.0, 10.0]))
    mask = numpy.array([False, False])
    fpeaks = mpu.filterPeaks(peaks, mask)

    assert(fpeaks.shape[0] == 0)


def test_load_mappings_1():
    map_test_file = storm_analysis.getPathOutputTest("map.map")

    max_ch = 4
    mappings = {}
    for i in range(1,max_ch):
        j = i
        mappings[str(i) + "_0_x"] = numpy.arange(j,j+2.5,1.0)
        j += 0.1
        mappings[str(i) + "_0_y"] = numpy.arange(j,j+2.5,1.0)
        j += 0.1
        mappings["0_" + str(i) + "_x"] = numpy.arange(j,j+2.5,1.0)
        j += 0.1
        mappings["0_" + str(i) + "_y"] = numpy.arange(j,j+2.5,1.0)

    max_ch -= 1

    with open(map_test_file, 'wb') as fp:
        pickle.dump(mappings, fp)

    mappings = {}
    [xt_0toN, yt_0toN, xt_Nto0, yt_Nto0] = mpUtilC.loadMappings(map_test_file, 0)
    assert(xt_0toN[0,0] == 0.0)
    assert(yt_0toN[0,0] == 0.0)
    assert(xt_Nto0[0,0] == 0.0)
    assert(yt_Nto0[0,0] == 0.0)

    assert(abs(xt_0toN[max_ch,2]-5.2) < 1.0e-6)
    assert(abs(yt_0toN[max_ch,2]-5.3) < 1.0e-6)
    assert(abs(xt_Nto0[max_ch,2]-5.0) < 1.0e-6)
    assert(abs(yt_Nto0[max_ch,2]-5.1) < 1.0e-6)

    
def test_mark_close_peaks_1():
    mpu = createMPU(1)
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 10.0, 10.0]),
                                numpy.array([10.0, 10.0, 20.0]))

    # Mark as converged.
    i_s = utilC.getStatusIndex()
    peaks[:,i_s] = 1.0

    # Set heights.
    i_h = utilC.getHeightIndex()
    peaks[:,i_h] = 10.0
    peaks[1,i_h] = 5.0
    
    [m_peaks, m_mask] = mpu.markClosePeaks(peaks)
    assert(m_peaks[0,i_s] == 0.0)
    assert(m_peaks[2,i_s] == 1.0)
    assert(m_mask[0] == 1)
    assert(m_mask[1] == 0)
    assert(m_mask[2] == 1)

        
def test_merge_new_peaks_1():
    # no overlap test.
    mpu = createMPU(1)
    peaks1 = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                 numpy.array([10.0, 10.0]))
    peaks2 = mpu.testCreatePeaks(numpy.array([20.0]),
                                 numpy.array([20.0]))
    m_peaks = mpu.mergeNewPeaks(peaks1, peaks2)
    assert(m_peaks.shape[0] == 3)

def test_merge_new_peaks_2():
    # complete overlap test.
    mpu = createMPU(1)
    peaks1 = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                 numpy.array([10.0, 10.0]))
    peaks2 = mpu.testCreatePeaks(numpy.array([10.0]),
                                 numpy.array([10.0]))
    i_s = utilC.getStatusIndex()
    peaks1[:,i_s] = 1.0
    m_peaks = mpu.mergeNewPeaks(peaks1, peaks2)
    assert(m_peaks.shape[0] == 2)

def test_merge_new_peaks_3():
    # partial overlap test.
    mpu = createMPU(1)
    peaks1 = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                 numpy.array([10.0, 10.0]))
    peaks2 = mpu.testCreatePeaks(numpy.array([12.0]),
                                 numpy.array([10.0]))
    i_s = utilC.getStatusIndex()
    peaks1[:,i_s] = 1.0
    m_peaks = mpu.mergeNewPeaks(peaks1, peaks2)
    assert(m_peaks.shape[0] == 3)
    assert(m_peaks[0,i_s] == 0.0)
    assert(m_peaks[1,i_s] == 1.0)


def test_remove_close_peaks_1():
    mpu = createMPU(1)
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 10.0, 10.0, 10.0]),
                                numpy.array([10.0, 10.0, 12.0, 20.0]))

    # Mark as converged.
    i_s = utilC.getStatusIndex()
    peaks[:,i_s] = 1.0

    # Set heights.
    i_h = utilC.getHeightIndex()
    peaks[:,i_h] = 10.0
    peaks[1,i_h] = 5.0
    
    peaks = mpu.removeClosePeaks(peaks)
    assert(peaks.shape[0] == 3)
    assert(peaks[0,i_s] == 0.0)
    assert(peaks[1,i_s] == 0.0)
    assert(peaks[2,i_s] == 1.0)

def test_remove_close_peaks_2():
    mpu = createMPU(2)
    mpu.setTransforms(*createTransforms(2))
    
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 10.0, 10.0, 10.0]),
                                numpy.array([10.0, 10.0, 12.0, 20.0]))
    peaks = mpu.splitPeaks(peaks)

    # Mark as converged.
    i_s = utilC.getStatusIndex()
    peaks[:,i_s] = 1.0

    # Set heights.
    i_h = utilC.getHeightIndex()
    peaks[:,i_h] = 10.0
    peaks[1,i_h] = 5.0
    
    peaks = mpu.removeClosePeaks(peaks)
    
    assert(peaks.shape[0] == 6)
    assert(peaks[0,i_s] == 0.0)
    assert(peaks[1,i_s] == 0.0)
    assert(peaks[2,i_s] == 1.0)
    assert(peaks[3,i_s] == 0.0)
    assert(peaks[4,i_s] == 0.0)
    assert(peaks[5,i_s] == 1.0)

    
def test_split_peaks_1():
    mpu = createMPU(1)
    mpu.setTransforms(*createTransforms(1))
    
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 20.0]),
                                numpy.array([10.0, 10.0]))
    s_peaks = mpu.splitPeaks(peaks)

    assert(s_peaks.shape == peaks.shape)
    peaks = peaks.flatten()
    s_peaks = s_peaks.flatten()
    for i in range(peaks.size):
        assert(peaks[i] == s_peaks[i])

def test_split_peaks_2():
    mpu = createMPU(2)
    mpu.setTransforms(*createTransforms(2))
    
    peaks = mpu.testCreatePeaks(numpy.array([10.0, 20.0, 30.0, 40.0]),
                                numpy.array([10.0, 10.0, 10.0, 10.0]))
    s_peaks = mpu.splitPeaks(peaks)

    assert(2 * peaks.shape[0] == s_peaks.shape[0])
    
    for i in range(peaks.shape[0]):
        assert(peaks[i,0] == s_peaks[i+peaks.shape[0],0])

    
if (__name__ == "__main__"):
    numpy.set_printoptions(precision=1, suppress=True)
    
    test_bad_peak_mask_1()
    test_bad_peak_mask_2()
    test_filter_1()
    test_filter_2()
    test_filter_3()
    test_filter_4()
    test_load_mappings_1()
    test_mark_close_peaks_1()
    test_merge_new_peaks_1()
    test_merge_new_peaks_2()
    test_merge_new_peaks_3()
    test_remove_close_peaks_1()
    test_remove_close_peaks_2()
    test_split_peaks_1()
    test_split_peaks_2()
