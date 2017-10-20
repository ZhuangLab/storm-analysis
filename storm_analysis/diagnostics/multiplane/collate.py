#!/usr/bin/env python
"""
Collate analysis results for Multiplane testing.

Hazen 10/17
"""
import glob
import math
import numpy

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.ia_utilities_c as utilC
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

import storm_analysis.sa_utilities.finding_fitting_error as ffe
import storm_analysis.sa_utilities.recall_fraction as rfrac

import settings

dirs = sorted(glob.glob("test*"))

if(len(dirs) == 0):
    print("No test directories found.")
    exit()

all_dx = []
all_dy = []
all_dz = []
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

    # Adjust channel 0 z values based on offset.
    i3_locs = readinsight3.loadI3File(a_dir + "/test_mlist.bin")
    zf = i3_locs['z'] + settings.z_planes[0]
    i3dtype.posSet(i3_locs, 'z', zf)
    with writeinsight3.I3Writer(a_dir + "/test_c1_mlist.bin") as i3w:
        i3w.addMolecules(i3_locs)
    
    # Load localizations.
    truth_i3 = readinsight3.I3Reader(a_dir + "/test_c1_olist.bin")
    measured_i3 = readinsight3.I3Reader(a_dir + "/test_c1_mlist.bin")
    total_locs += measured_i3.getNumberMolecules()

    # Calculate fractional recall.
    [partial, total] = rfrac.recallFraction(truth_i3, measured_i3, settings.tolerance)
    recall += partial
    recall_total += total

    # Calculate noise fraction.
    [partial, total] = rfrac.noiseFraction(truth_i3, measured_i3, settings.tolerance)
    noise += partial
    noise_total += total

    # Calculate fitting error in XYZ.
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
        all_dz.append([numpy.std(dz), math.sqrt(numpy.mean(dz*dz))])
    else:
        all_dx.append([0,0])
        all_dy.append([0,0])
        all_dz.append([0,0])


print()
print("Analysis Summary:")
print("Processed {0:0d} localizations in {1:.2f} seconds, {2:.2f}/sec".format(total_locs, total_time, float(total_locs)/float(total_time)))
print("Recall {0:.5f}".format(float(recall)/float(recall_total)))
print("Noise {0:.5f}".format(float(noise)/float(noise_total)))
print("XYZ Precision (nm):")
for i, a_dir in enumerate(dirs):
    print(a_dir + "\t{0:.2f}\t{1:.2f}\t{2:.2f}".format(all_dx[i][0], all_dy[i][0], all_dz[i][0]))
print("")
print("XYZ RMS Accuracy (nm):")
for i, a_dir in enumerate(dirs):
    print(a_dir + "\t{0:.2f}\t{1:.2f}\t{2:.2f}".format(all_dx[i][1], all_dy[i][1], all_dz[i][1]))
print("")
