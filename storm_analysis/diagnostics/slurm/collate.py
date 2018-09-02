#!/usr/bin/env python
"""
Collate analysis results for SLURM testing.

Hazen 09/18
"""
import glob
import os

import storm_analysis.slurm.check_analysis as checkAnalysis
import storm_analysis.slurm.merge_analysis as mergeAnalysis

import storm_analysis.diagnostics.slurm.settings as settings


def collate():

    # Remove stale results, if any.
    m_name = os.path.join(settings.wdir, "merged.hdf5")
    if os.path.exists(m_name):
        os.remove(m_name)
            
    checkAnalysis.checkAnalysis(settings.wdir)
    mergeAnalysis.mergeAnalysis(settings.wdir, m_name)
    

if (__name__ == "__main__"):
    collate()
