#!/usr/bin/env python
"""
Collate analysis results.

Hazen 09/17
"""
import glob
import math
import numpy

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC

import storm_analysis.sa_utilities.finding_fitting_error as ffe
import storm_analysis.sa_utilities.recall_fraction as rfrac

import settings

dirs = sorted(glob.glob("test*"))

if(len(dirs) == 0):
    print("No test directories found.")
    exit()

all_dx = []
all_dy = []
all_wx = []
all_wy = []
noise = 0
noise_total = 0
recall = 0
recall_total = 0
total_locs = 0
total_time = 0.0
for a_dir in dirs:
    print("Processing", a_dir)

    # Load timing information.
    with open(a_dir + "/timing.txt") as fp:
        total_time += float(fp.readline())

    # Load localizations.
    truth_i3 = readinsight3.I3Reader(a_dir + "/test_olist.bin")
    measured_i3 = readinsight3.I3Reader(a_dir + "/test_mlist.bin")
    total_locs += measured_i3.getNumberMolecules()
    
    # Calculate fractional recall.
    [partial, total] = rfrac.recallFraction(truth_i3, measured_i3, settings.tolerance)
    recall += partial
    recall_total += total

    # Calculate noise fraction.
    [partial, total] = rfrac.noiseFraction(truth_i3, measured_i3, settings.tolerance)
    noise += partial
    noise_total += total

    # Calculate error in fitting width.
    for i in range(truth_i3.getNumberFrames()):    
        t_locs = truth_i3.getMoleculesInFrame(i+1)
        m_locs = measured_i3.getMoleculesInFrame(i+1)

        # Widths for truth localizations.
        ax = t_locs['ax']
        ww = t_locs['w']
        t_wx = 0.5*numpy.sqrt(ww*ww/ax)/settings.pixel_size
        t_wy = 0.5*numpy.sqrt(ww*ww*ax)/settings.pixel_size
    
        # Widths for found localizations.
        ax = m_locs['ax']
        ww = m_locs['w']
        m_wx = 0.5*numpy.sqrt(ww*ww/ax)/settings.pixel_size
        m_wy = 0.5*numpy.sqrt(ww*ww*ax)/settings.pixel_size
        
        p_index = iaUtilsC.peakToPeakDistAndIndex(m_locs['xc'], m_locs['yc'],
                                                  t_locs['xc'], t_locs['yc'],
                                                  max_distance = settings.tolerance)[1]

        p_size = numpy.count_nonzero(p_index > -1)
        d_wx = numpy.zeros(p_size)
        d_wy = numpy.zeros(p_size)
        k = 0
        for j in range(m_locs.size):
            if(p_index[j] < 0):
                continue
            
            d_wx[k] = m_wx[j] - t_wx[p_index[j]]
            d_wy[k] = m_wy[j] - t_wy[p_index[j]]
            k += 1
        
        all_wx.append(d_wx)
        all_wy.append(d_wy)

    # Calculate fitting error in XY.
    max_distance = None
    if True:
        max_distance = 2.0 * settings.pixel_size
        print("Using max_distance", max_distance, "nm for error calcuations.")
        
    [dx, dy, dz] = ffe.findingFittingError(truth_i3,
                                           measured_i3,
                                           pixel_size = settings.pixel_size,
                                           max_distance = max_distance)
    if dx is not None:
        all_dx.append([numpy.std(dx), math.sqrt(numpy.mean(dx*dx))])
        all_dy.append([numpy.std(dy), math.sqrt(numpy.mean(dy*dy))])
    else:
        all_dx.append([0,0])
        all_dy.append([0,0])


print()
print("Analysis Summary:")
print("Processed {0:0d} localizations in {1:.2f} seconds, {2:.2f}/sec".format(total_locs, total_time, float(total_locs)/float(total_time)))
print("Recall {0:.5f}".format(float(recall)/float(recall_total)))
print("Noise {0:.5f}".format(float(noise)/float(noise_total)))
print("XY Precision (nm):")
for i, a_dir in enumerate(dirs):
    print(a_dir + "\t{0:.2f}\t{1:.2f}".format(all_dx[i][0], all_dy[i][0]))
print("")
print("XY RMS Accuracy (nm):")
for i, a_dir in enumerate(dirs):
    print(a_dir + "\t{0:.2f}\t{1:.2f}".format(all_dx[i][1], all_dy[i][1]))
print("")
print("XY Width Error, Mean difference with truth, Standard deviation (pixels):")
for i, a_dir in enumerate(dirs):
    print(a_dir + "\t{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}".format(numpy.mean(all_wx[i]),
                                                                numpy.std(all_wx[i]),
                                                                numpy.mean(all_wy[i]),
                                                                numpy.std(all_wy[i])))

