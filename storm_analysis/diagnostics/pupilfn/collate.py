#!/usr/bin/env python
"""
Collate analysis results for Pupilfn testing.

Hazen 10/17
"""
import glob

import storm_analysis.diagnostics.collate as collateResults

import storm_analysis.diagnostics.pupilfn.settings as settings


def collate():
    dirs = sorted(glob.glob("test*"))

    if(len(dirs) == 0):
        print("No test directories found.")
        return

    collateResults.collateSpliner(dirs, settings)


def collateRQE():
    dirs = sorted(glob.glob("test*"))

    if(len(dirs) == 0):
        print("No test directories found.")
        return

    collateResults.collateRQE(dirs, settings)
    
    
if (__name__ == "__main__"):
    collate()
    
