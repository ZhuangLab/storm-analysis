#!/usr/bin/env python
"""
Collate analysis results for Multicolor testing.

Hazen 01/18
"""
import glob
import inspect
import numpy

import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.diagnostics.multicolor.settings as settings

def collate():

    dirs = sorted(glob.glob("test*"))

    if(len(dirs) == 0):
        print("No test directories found.")
        exit()

    # Load reference localizations.
    ref = saH5Py.loadLocalizations("sim_input_c1_grid_list.hdf5", fields = ["color", "x", "y"])
    
    for a_dir in dirs:
    
        # Check color correspondence.
        #
        # Note: Currently only works for gridded localizations.
        #
        exp = saH5Py.loadTracks(a_dir + "/test.hdf5", fields = ["x", "y", "km_color"])

        # Identify corresponding peaks with 2 pixel maximum radius.
        p_index = iaUtilsC.peakToPeakDistAndIndex(exp["x"], exp["y"], ref["x"], ref["y"],
                                                  max_distance = 2.0)[1]

        # Create array for checking category correspondence.
        c_size = numpy.count_nonzero(p_index > -1)
        categories = numpy.zeros((c_size, 2), dtype = numpy.int32)
        i = 0
        for j in range(p_index.size):
            if (p_index[j] < 0):
                continue
            categories[i,0] = int(ref["color"][p_index[j]])
            categories[i,1] = int(exp["km_color"][j])
            i += 1

        # Measure fraction of experimental localizations assigned correctly. This
        # is complicated by the k-mean category being arbitrary.
        for i in range(8):
            ref_mask = (categories[:,0] == i)
            if(numpy.count_nonzero(ref_mask) > 0):
                exp_cat = int(numpy.median(categories[ref_mask,1]))
                exp_mask = (categories[:,1] == exp_cat)

                total = numpy.count_nonzero(exp_mask)
                matched = numpy.count_nonzero(numpy.logical_and(ref_mask, exp_mask))
                mismatched = total - matched
                print("Color {0:0d}, matching fraction is {1:.2f}, unmatched fraction is {2:.2f}, total {3:0d}".format(i, matched/total, mismatched/total, total))
        
if (__name__ == "__main__"):
    collate()
