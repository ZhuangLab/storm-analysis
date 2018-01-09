#!/usr/bin/env python
"""
Collate analysis results for Multiplane testing.

Hazen 10/17
"""
import glob

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.diagnostics.collate as collate

import settings

dirs = sorted(glob.glob("test*"))

if(len(dirs) == 0):
    print("No test directories found.")
    exit()

# Adjust z positions in the channel 0 reference.
for a_dir in dirs:
    with saH5Py.SAH5Py(a_dir + "/test_c1_ref.hdf5") as h5_in:
        with saH5Py.SAH5Py(a_dir + "/test_ref.hdf5", is_existing = False, overwrite = True) as h5_out:
            h5_out.setMovieInformation(*h5_in.getMovieInformation())

            for fnum, locs in h5_in.localizationsIterator():
                locs["z"] -= 1.0e-3 * settings.z_planes[0]
                h5_out.addLocalizations(locs, fnum)

# Collate results.
collate.collateSpliner(dirs, settings)
