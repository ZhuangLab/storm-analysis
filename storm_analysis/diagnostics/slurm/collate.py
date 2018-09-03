#!/usr/bin/env python
"""
Collate analysis results for SLURM testing.

Hazen 09/18
"""
import glob
import os

import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.slurm.check_analysis as checkAnalysis
import storm_analysis.slurm.merge_analysis as mergeAnalysis

import storm_analysis.diagnostics.slurm.settings as settings


def collate():

    # Remove stale results, if any.
    m_name = os.path.join(settings.wdir, "merged.hdf5")
    if os.path.exists(m_name):
        os.remove(m_name)

    # Run SLURM functions.
    checkAnalysis.checkAnalysis(settings.wdir)
    mergeAnalysis.mergeAnalysis(settings.wdir, m_name)

    # Check results.
    with saH5Py.SAH5Py(m_name) as h5:
        n_locs_slurm = h5.getNLocalizations()

    n_locs = 0
    h5_names = sorted(glob.glob(os.path.join(settings.wdir, "p_*.hdf5")))
    for h5_name in h5_names:
        with saH5Py.SAH5Py(h5_name) as h5:
            n_locs += h5.getNLocalizations()        

    if (n_locs_slurm == n_locs):
        print("All good.")
    else:
        print("Localization counts do not match", n_locs, n_locs_slurm)
    

if (__name__ == "__main__"):
    collate()
