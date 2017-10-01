#!/usr/bin/env python
"""
Collate analysis results for Spliner testing.

Hazen 09/17
"""
import glob
import numpy

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.ia_utilities_c as utilC

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
total_time = 0.0
for a_dir in dirs:
    print("Processing", a_dir)

    # Load timing information.
    with open(a_dir + "/timing.txt") as fp:
        total_time += float(fp.readline())

    # Load localizations.
    truth_i3 = readinsight3.I3Reader(a_dir + "/test_olist.bin")
    measured_i3 = readinsight3.I3Reader(a_dir + "/test_mlist.bin")    
    
    # Calculate fractional recall.
    [partial, total] = rfrac.recallFraction(truth_i3, measured_i3, settings.tolerance)
    recall += partial
    recall_total += total

    # Calculate noise fraction.
    [partial, total] = rfrac.noiseFraction(truth_i3, measured_i3, settings.tolerance)
    noise += partial
    noise_total += total

    # Calculate fitting error in XYZ.
    [dx, dy, dz] = ffe.findingFittingError(truth_i3, measured_i3, pixel_size = settings.pixel_size)
    if dx is not None:
        all_dx.append(numpy.std(dx))
        all_dy.append(numpy.std(dy))
        all_dz.append(numpy.std(dz))
    else:
        all_dx.append(0)
        all_dy.append(0)
        all_dz.append(0)


print()
print("Analysis Summary:")
print("Total analysis time {0:.2f} seconds".format(total_time))
print("Recall {0:.5f}".format(float(recall)/float(recall_total)))
print("Noise {0:.5f}".format(float(noise)/float(noise_total)))
print("XYZ Error (nm):")
for i, a_dir in enumerate(dirs):
    print(a_dir + "\t{0:.2f}\t{1:.2f}\t{2:.2f}".format(all_dx[i], all_dy[i], all_dz[i]))

