#!/usr/bin/env python
"""
Collate analysis results.

Hazen 09/17
"""
import glob

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.diagnostics.collate as collateResults

import storm_analysis.diagnostics.multiplane_dao.settings as settings


def collate():
    dirs = sorted(glob.glob("test*"))

    if(len(dirs) == 0):
        print("No test directories found.")
        return

    # Adjust z positions in the channel 0 reference.
    for a_dir in dirs:
        with saH5Py.SAH5Py(a_dir + "/test_c1_ref.hdf5") as h5_in:
            with saH5Py.SAH5Py(a_dir + "/test_ref.hdf5", is_existing = False, overwrite = True) as h5_out:
                h5_out.setMovieInformation(*h5_in.getMovieInformation())

                for fnum, locs in h5_in.localizationsIterator():
                    locs["z"] -= settings.z_planes[0]
                    h5_out.addLocalizations(locs, fnum)
                
    collateResults.collateDAO(dirs, settings, calc_width_error = False)


if (__name__ == "__main__"):
    collate()
    
